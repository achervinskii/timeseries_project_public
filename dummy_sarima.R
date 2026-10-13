## =====================================================================
#  SARIMA Rolling-Window Forecast + Diebold Mariano Test
## =====================================================================

library(ggplot2)
library(lubridate)
library(dplyr)
library(tidyr)
library(purrr)
library(tsibble)
library(fable)
library(feasts)
library(fabletools)
library(forecast)
library(MCS)

## ----- Loading Data -----
load("citibike.RData")

df <- data_citibike
df$datetime <- ymd_h(paste(df$year, df$month, df$day, df$hour))
df$weekday <- wday(df$datetime, label = TRUE)
df$is_weekend <- factor(ifelse(df$weekday %in% c("Sat", "Sun"), "Weekend", "Workday"))
names(df) <- gsub(" ", "", names(df))

## ----- Training sample: April until July 2023 -----
train <- df %>% filter(datetime < ymd("2023-08-01"))
test <- df %>% filter(datetime >= ymd("2023-08-01"))

## ----- Build factors -----
df <- df %>%
  mutate(
    hour_f  = factor(hour, levels = 0:23),
    month_f = factor(month, levels = 1:12)
  )

train_reg <- train %>%
  mutate(
    hour_f  = factor(hour, levels = 0:23),
    month_f = factor(month, levels = 1:12)
  )

train_reg_tsbl <- train_reg %>%
  select(datetime, area6, is_weekend, hour_f, month_f) %>%
  as_tsibble(index = datetime)

## ----- OLS weekend/workday x hour interaction (exploratory fit on full training data) -----
fit_reg <- train_reg_tsbl %>%
  model(reg = TSLM(area6 ~ 0 + is_weekend:hour_f + month_f))

report(fit_reg)
tidy(fit_reg)

## ----- Extract epsilon_t for the SARIMA step -----
resid_tsbl <- residuals(fit_reg) %>%
  select(datetime, .resid) %>%
  rename(epsilon = .resid) %>%
  as_tsibble(index = datetime)

plt <- resid_tsbl %>% autoplot(epsilon) +
  labs(title = "Residuals from the weekend/workday x hour regression", y = "epsilon")
ggsave(
  plt,
  filename = file.path("outputs", "residuals_from_dummy_regression.png"),
  width = 10, height = 6
)

## ----- Correlogram of epsilon_t -----
png("outputs/partial_corr_eps.png", width = 1200, height = 800, res = 150)
plt <- resid_tsbl %>% gg_tsdisplay(epsilon, plot_type = "partial", lag_max = 100)

ggsave(
  plt,
  filename = file.path("outputs", "partial_corr_eps.png"),
  width = 10, height = 6
)

## ----- Fit SARIMA on epsilon_t -----
# d = 0, D = 0 (no differencing required - trend stationary)
fit_sarima <- resid_tsbl %>%
  model(
    sarima_1    = ARIMA(epsilon ~ pdq(4, 0, 0) + PDQ(0, 0, 0, period = 24)),
    sarima_2    = ARIMA(epsilon ~ pdq(4, 0, 0) + PDQ(1, 0, 0, period = 24)),
    # sarima_auto = ARIMA(epsilon, stepwise = FALSE, approximation = FALSE)
  )

glance(fit_sarima) %>% select(.model, AIC, AICc, BIC)

fit_sarima %>%
  select(sarima_1) %>%
  report()
fit_sarima %>%
  select(sarima_2) %>%
  report()
# fit_sarima %>%
#   select(sarima_auto) %>%
#   report()


png("outputs/sarima_fit_residuals.png", width = 1200, height = 800, res = 150)
fit_sarima %>%
  select(sarima_1) %>%
  gg_tsresiduals()
dev.off()

png("outputs/sarima_2_fit_residuals.png", width = 1200, height = 800, res = 150)
fit_sarima %>%
  select(sarima_2) %>%
  gg_tsresiduals()
dev.off()

# fit_sarima %>%
#   select(sarima_auto) %>%
#   gg_tsresiduals()
#
## ----- Ljung-Box test -----
dofs <- c(sarima_1 = 3, sarima_2 = 5)
# dofs <- c(sarima_1 = 3, sarima_2 = 5, sarima_auto = 5)

purrr::imap_dfr(dofs, \(d, m)
augment(fit_sarima) %>%
  filter(.model == m) %>%
  features(.innov, ljung_box, lag = 48, dof = d))

## =====================================================================
## Rolling-window 24h-ahead SARIMA forecasts
## =====================================================================

sarima_spec <- ARIMA(epsilon ~ pdq(4, 0, 0) + PDQ(1, 0, 0, period = 24))

origins <- test %>%
  distinct(date = as_date(datetime)) %>%
  pull(date) %>%
  as_datetime() # one origin per calendar day, 00:00

window_length <- nrow(train) # rolling window length = size of original training sample (hours)

roll_results <- vector("list", length(origins))

for (i in seq_along(origins)) {
  origin <- origins[i]

  # --- rolling window: window_length hours immediately before `origin` ---
  window_data <- df %>%
    filter(datetime >= origin - hours(window_length), datetime < origin)

  window_tsbl <- window_data %>%
    select(datetime, area6, is_weekend, hour_f, month_f) %>%
    as_tsibble(index = datetime)

  # --- regression step ---
  # Fallback: only the very first origin (Aug 1) has a training window with
  # zero August observations, so month_f can't be estimated there. We drop
  # month_f for that one origin only; every other origin uses the full spec.
  fit_reg_i <- tryCatch(
    window_tsbl %>% model(reg = TSLM(area6 ~ 0 + is_weekend:hour_f + month_f)),
    warning = function(w) window_tsbl %>% model(reg = TSLM(area6 ~ 0 + is_weekend:hour_f)),
    error = function(e) window_tsbl %>% model(reg = TSLM(area6 ~ 0 + is_weekend:hour_f))
  )

  resid_i <- residuals(fit_reg_i) %>%
    select(datetime, .resid) %>%
    rename(epsilon = .resid) %>%
    as_tsibble(index = datetime)

  # --- SARIMA step on epsilon (fixed order, no re-search) ---
  fit_sarima_i <- resid_i %>% model(sarima = sarima_spec)
  fc_epsilon_i <- fit_sarima_i %>% forecast(h = 24)

  # --- next 24h regressors ---
  newdata_i <- df %>%
    filter(datetime >= origin, datetime < origin + hours(24)) %>%
    select(datetime, area6, is_weekend, hour_f, month_f) %>%
    as_tsibble(index = datetime)

  fc_mean_i <- fit_reg_i %>% forecast(new_data = newdata_i)

  roll_results[[i]] <- tibble(
    datetime = newdata_i$datetime,
    horizon  = 1:24,
    actual   = newdata_i$area6,
    forecast = fc_mean_i$.mean + fc_epsilon_i$.mean
  )
}

roll_fc <- bind_rows(roll_results)

# Clip negative forecasts at zero FOR PLOTTING ONLY — bike demand can't be
# negative, but we keep the raw 'forecast' column intact for MAE / DM test,
# so reported accuracy isn't altered by the clipping.
roll_fc <- roll_fc %>% mutate(forecast_plot = pmax(forecast, 0))

plt <- ggplot(roll_fc, aes(x = datetime)) +
  geom_line(aes(y = actual, colour = "Actual")) +
  geom_line(aes(y = forecast_plot, colour = "SARIMA rolling forecast")) +
  labs(
    title = "Rolling 24h-ahead forecasts vs actual demand (area 6, August 2023)",
    x = NULL, y = "Demand", colour = NULL
  ) +
  theme_minimal()

ggsave(
  plt,
  filename = file.path("outputs", "rolling_vs_actual.png"),
  width = 10, height = 6
)
## ----- MAE -----
mae_overall <- roll_fc %>%
  summarise(MAE = mean(abs(actual - forecast), na.rm = TRUE))

mae_by_horizon_SARIMA <- roll_fc %>%
  group_by(horizon) %>%
  summarise(MAE = mean(abs(actual - forecast), na.rm = TRUE))

print(mae_by_horizon_SARIMA, n = 24)

## ----- Horizon-averaged forecast vs actual (SARIMA) -----
horizon_avg <- roll_fc %>%
  group_by(horizon) %>%
  summarise(
    forecast_avg = mean(forecast, na.rm = TRUE),
    actual_avg   = mean(actual, na.rm = TRUE)
  )

results_sarima <- horizon_avg %>%
  rename(h = horizon, SARIMA = forecast_avg, actual = actual_avg)


png("outputs/sarima_forecast_vs_actual.png", width = 1200, height = 800, res = 150)
matplot(results_sarima$h, results_sarima[, c("SARIMA", "actual")],
  type = "l", lty = 1,
  col = c("magenta", "turquoise"), xlab = "Horizon h", ylab = "Demand"
)
legend("topleft", c("SARIMA Forecast", "Actual"), col = c("magenta", "turquoise"), lty = 1)
dev.off()

## =====================================================================
## Diebold-Mariano test: SARIMA vs ETS
## =====================================================================

## ----- Build SARIMA raw error matrix (31 origins x 24 horizons) -----
err_sarima <- roll_fc %>%
  mutate(
    origin_id = rep(1:(n() / 24), each = 24),
    error = actual - forecast
  ) %>%
  select(origin_id, horizon, error) %>%
  pivot_wider(id_cols = origin_id, names_from = horizon, values_from = error) %>%
  select(-origin_id) %>%
  as.matrix()

dim(err_sarima)

## ----- Load ETS raw error vector and reshape -----
load("outputs/err_ETS.RData")

err_ETS_mat <- matrix(err_ETS, nrow = 31, ncol = 24, byrow = TRUE)
dim(err_ETS_mat) # expect 31 x 24

# Save the properly shaped MATRIX (not the flat vector) for future reuse
save(err_ETS_mat, file = "outputs/err_ETS_matrix.RData")

## ----- DM test per horizon -----
dm_results <- map_dfr(1:24, function(h) {
  test <- dm.test(err_sarima[, h], err_ETS_mat[, h],
    alternative = "two.sided", h = h, power = 1,
    varestimator = "bartlett"
  )
  tibble(horizon = h, DM_stat = as.numeric(test$statistic), p_value = test$p.value)
})

print(dm_results, n = 24)

## =====================================================================
## Combined MAE-by-horizon plot with DM significance markers
## =====================================================================

mae_by_horizon_ETS <- tibble(
  horizon = 1:24,
  MAE = colMeans(abs(err_ETS_mat), na.rm = TRUE)
)

mae_compare <- bind_rows(
  mae_by_horizon_SARIMA %>% mutate(model = "SARIMA"),
  mae_by_horizon_ETS %>% mutate(model = "ETS")
)

sig_horizons <- dm_results %>%
  mutate(
    significant = p_value < 0.05,
    winner = case_when(
      !significant ~ "None",
      DM_stat < 0 ~ "SARIMA",
      DM_stat > 0 ~ "ETS"
    )
  ) %>%
  select(horizon, significant, winner)

mae_compare <- mae_compare %>% left_join(sig_horizons, by = "horizon")

marker_y <- mae_compare %>%
  group_by(horizon) %>%
  summarise(y = max(MAE) * 1.08, winner = first(winner), significant = first(significant))

plt <- ggplot(mae_compare, aes(x = horizon, y = MAE, colour = model)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  geom_point(
    data = marker_y %>% filter(significant),
    aes(x = horizon, y = y, colour = winner),
    shape = 8, size = 2.5, show.legend = FALSE
  ) +
  scale_colour_manual(values = c("SARIMA" = "#1b9e77", "ETS" = "#d95f02")) +
  labs(
    title = "MAE by forecast horizon: SARIMA vs ETS (area 6, August 2023)",
    subtitle = "Asterisk marks horizons with a significant DM test (5% level)",
    x = "Horizon (hours ahead)", y = "MAE", colour = NULL
  ) +
  theme_minimal(base_size = 12)

ggsave(
  plt,
  filename = file.path("outputs", "MAE_sarima_vs_ets.png"),
  width = 10, height = 6
)


## ---- MCS Test ---- 
# MCS test for multiple horizons (1-24) comparing SARIMA and ETS
loss_matrix <- cbind(
  SARIMA = abs(as.vector(err_sarima)),
  ETS    = abs(as.vector(err_ETS_mat))
)

mcs_result <- MCSprocedure(Loss = loss_matrix, 
                           alpha = 0.05, 
                           B = 5000,
                           statistic = "Tmax", 
                           verbose = TRUE)


## ---- MCS TEST WITH 3 MODELS ----
sarima_spec_1 <- ARIMA(epsilon ~ pdq(4, 0, 0) + PDQ(0, 0, 0, period = 24))

roll_results_1 <- vector("list", length(origins))

for (i in seq_along(origins)) {
  
  origin <- origins[i]
  
  window_data <- df %>%
    filter(datetime >= origin - hours(window_length), datetime < origin)
  
  window_tsbl <- window_data %>%
    select(datetime, area6, is_weekend, hour_f, month_f) %>%
    as_tsibble(index = datetime)
  
  fit_reg_i <- tryCatch(
    window_tsbl %>% model(reg = TSLM(area6 ~ 0 + is_weekend:hour_f + month_f)),
    warning = function(w) window_tsbl %>% model(reg = TSLM(area6 ~ 0 + is_weekend:hour_f)),
    error   = function(e) window_tsbl %>% model(reg = TSLM(area6 ~ 0 + is_weekend:hour_f))
  )
  
  resid_i <- residuals(fit_reg_i) %>%
    select(datetime, .resid) %>%
    rename(epsilon = .resid) %>%
    as_tsibble(index = datetime)
  
  fit_sarima_i <- resid_i %>% model(sarima = sarima_spec_1)   # <- sarima_1 spec here
  fc_epsilon_i <- fit_sarima_i %>% forecast(h = 24)
  
  newdata_i <- df %>%
    filter(datetime >= origin, datetime < origin + hours(24)) %>%
    select(datetime, area6, is_weekend, hour_f, month_f) %>%
    as_tsibble(index = datetime)
  
  fc_mean_i <- fit_reg_i %>% forecast(new_data = newdata_i)
  
  roll_results_1[[i]] <- tibble(
    datetime = newdata_i$datetime,
    horizon  = 1:24,
    actual   = newdata_i$area6,
    forecast = fc_mean_i$.mean + fc_epsilon_i$.mean
  )
}

roll_fc_1 <- bind_rows(roll_results_1)

err_sarima_1 <- roll_fc_1 %>%
  mutate(origin_id = rep(1:(n() / 24), each = 24),
         error = actual - forecast) %>%
  select(origin_id, horizon, error) %>%
  pivot_wider(id_cols = origin_id, names_from = horizon, values_from = error) %>%
  select(-origin_id) %>%
  as.matrix()

dim(err_sarima_1)   # should be 31 x 24, matching err_sarima and err_ETS_mat

loss_matrix_3 <- cbind(
  SARIMA_2 = abs(as.vector(err_sarima)),      # your original ARIMA(4,0,0)(1,0,0)[24]
  SARIMA_1 = abs(as.vector(err_sarima_1)),    # ARIMA(4,0,0), no seasonal term
  ETS      = abs(as.vector(err_ETS_mat))
)

mcs_result_3 <- MCSprocedure(Loss = loss_matrix_3, alpha = 0.05, B = 5000,
                             statistic = "Tmax", verbose = TRUE)
