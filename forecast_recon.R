library(ggplot2)
library(lubridate)
library(tidyr)
library(dplyr)
library(tsibble)
library(fable)
library(feasts)
library(fabletools)
library(forecast)
library(FoReco)

## ----- Loading Data -----
load("citibike.RData")

df <- data_citibike
names(df) <- gsub(" ", "", names(df))

# colnames of areas
areas <- paste0("area", 1:9)
# compute all the factors
df_selected <- df |>
  mutate(
    datetime = ymd_h(paste(year, month, day, hour)),
    weekday = wday(datetime, label = TRUE),
    is_weekend_f = factor(ifelse(weekday %in% c("Sat", "Sun"),
      "Weekend", "Workday"
    )),
    hour_f = factor(hour, levels = 0:23),
    hour_weekend = interaction(
      is_weekend_f,
      hour_f,
      drop = TRUE,
      sep = ":"
    ),
    month_f = factor(month, levels = 1:12),
    area7 = area1 + area2 + area3 + area4,
    area8 = area5 + area6,
    area9 = area7 + area8
  ) |>
  select(datetime, month_f, hour_f, hour_weekend, all_of(areas))

head(df_selected)

# compute mean hourly demand for every area, rounded to the integer
mean_demands <- sapply(areas,
  FUN = function(i) round(mean(df_selected[[i]]), 0)
)
mean_demands

## ----- Training sample: April until July 2023 -----
train <- df_selected |> filter(datetime < ymd("2023-08-01"))
test <- df_selected |> filter(datetime >= ymd("2023-08-01"))

## ----- Rolling-window 24h-ahead SARIMA forecasts -----
sarima_spec <- ARIMA(epsilon ~ pdq(4, 0, 0) + PDQ(1, 0, 0, period = 24))

# origin is the start of a new day
origins <- test |>
  distinct(date = as_date(datetime)) |>
  pull(date) |>
  as_datetime() # one origin per calendar day, 00:00

# rolling window = length of original training sample (hours)
window_length <- nrow(train)

# arrays to keep the forecasts
base_forecasts <- array(NA, dim = c(length(areas), length(origins), 24))
reconciled_forecasts <- array(NA, dim = c(length(areas), length(origins), 24))
bu_forecasts <- array(NA, dim = c(length(areas), length(origins), 24))
shr_forecasts <- array(NA, dim = c(length(areas), length(origins), 24))
realizations <- array(NA, dim = c(length(areas), length(origins), 24))

# reconciliation matrix
agg_mat <- matrix(
  c(
    1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 0, 0,
    0, 0, 0, 0, 1, 1
  ),
  nrow = 3,
  byrow = TRUE
)
options(warn = 0)
for (i in seq_along(origins)) {
  origin <- origins[i]

  # --- rolling window: window_length hours immediately before `origin` ---
  window_data <- df_selected |>
    filter(datetime >= origin - hours(window_length), datetime < origin)

  # residual matrix
  resid_matrix <- matrix(NA, nrow(window_data), length(areas))

  # forecast for every area
  for (a in seq_along(areas)) {
    window_tsbl <- window_data |>
      select(datetime,
        y = all_of(areas[a]), hour_weekend, month_f
      ) |>
      as_tsibble(index = datetime)

    # --- regression step ---
    fit_reg_i <- lm(y ~ hour_weekend + month_f, data = window_tsbl)
    resid_i <- residuals(fit_reg_i)

    # --- auto-SARIMA step on epsilon
    fit_sarima_i <- auto.arima(ts(resid_i, frequency = 24),
      max.p = 4,
      max.q = 2,
      max.P = 1,
      max.Q = 1,
      max.order = 6,
      seasonal = TRUE
    )
    fc_res_i <- forecast(fit_sarima_i, h = 24)$mean

    # --- in-sample fitted values for reconciliation ---
    reg_fitted <- fitted(fit_reg_i)
    sarima_fitted <- fitted(fit_sarima_i)
    final_fitted <- reg_fitted + sarima_fitted
    final_resid <- window_tsbl$y - final_fitted
    resid_matrix[, a] <- final_resid

    # --- next 24h regressors (known in advance) ---
    # because the month of forecasting was not yet fully seen in the training
    # decrment it by onne when the month is not available in the training data.

    forecast_month <- month(origin)
    available_months <- as.integer(levels(droplevels(window_tsbl$month_f)))

    forecast_month_used <- if (forecast_month %in% available_months) {
      forecast_month
    } else {
      previous_month <- if (forecast_month == 1) 12 else forecast_month - 1
      previous_month
    }

    newdata_i <- df_selected |>
      filter(datetime >= origin, datetime < origin + hours(24)) |>
      select(datetime, y = all_of(areas[a]), hour_weekend, month_f) |>
      as_tsibble(index = datetime)

    # update the month
    newdata_i$month_f[] <- forecast_month_used

    fc_reg_i <- predict(fit_reg_i, newdata = newdata_i)

    base_forecasts[a, i, 1:24] <- fc_reg_i + fc_res_i
    realizations[a, i, 1:24] <- newdata_i$y
    print(paste(
      "origin", i,
      "area", a
    ))
  }

  # go over every horizon, and reuse the residual matrix to reconcile
  for (h in 1:24) {
    # foreco package uses different ordering, with bottom series at the end
    order <- c(9, 7, 8, 1, 2, 3, 4, 5, 6)
    inv_order <- c(4, 5, 6, 7, 8, 9, 2, 3, 1)
    base_for <- matrix(base_forecasts[seq_along(areas), i, h], nrow = 1)
    # foreco uses different ordering
    rec_for <- csrec(base_for[, order, drop = FALSE],
      agg_mat = agg_mat,
      comb = "wls",
      res = resid_matrix[, order, drop = FALSE]
    )
    bu_for <- csbu(base_for[, 1:6, drop = FALSE],
      agg_mat = agg_mat
    )
    shr_for <- csrec(
      base = base_for[, order, drop = FALSE],
      agg_mat = agg_mat,
      comb = "shr",
      res = resid_matrix[, order, drop = FALSE]
    )
    reconciled_forecasts[seq_along(areas), i, h] <- rec_for[1, inv_order]
    bu_forecasts[seq_along(areas), i, h] <- bu_for[1, inv_order]
    shr_forecasts[seq_along(areas), i, h] <- shr_for[1, inv_order]
  }
}

# plot MAE for base_foraecasts, reconciled_forecasts, and bu_forecasts.
# different forecasts in different colors
# compute mean over origins
mae_base <- apply(abs(realizations - base_forecasts),
  MARGIN = c(1, 3), FUN = mean
)

mae_rec <- apply(abs(realizations - reconciled_forecasts),
  MARGIN = c(1, 3), FUN = mean
)

mae_bu <- apply(abs(realizations - bu_forecasts),
  MARGIN = c(1, 3), FUN = mean
)

mae_shr <- apply(abs(realizations - shr_forecasts),
  MARGIN = c(1, 3), FUN = mean
)

mae_df <- expand_grid(
  horizon = 1:24,
  area = areas
) |>
  mutate(
    Base = as.vector(mae_base),
    WLS = as.vector(mae_rec),
    Bottom_up = as.vector(mae_bu),
    Shrinkage = as.vector(mae_shr)
  ) |>
  pivot_longer(
    c(Base, WLS, Bottom_up, Shrinkage),
    names_to = "forecast",
    values_to = "MAE"
  )

plt <- ggplot(mae_df, aes(horizon, MAE, color = forecast)) +
  geom_line(linewidth = 0.4) +
  facet_wrap(~area, scales = "free_y") +
  labs(
    x = "Forecast horizon",
    y = "MAE",
    color = "Forecast"
  )

ggsave(
  plt,
  filename = file.path("outputs", "mae_reconciled.png"),
  width = 10, height = 6
)
