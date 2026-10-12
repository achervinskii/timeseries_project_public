# install.packages("bootUR")
# install.packages("tsibble")
# install.packages("fable")
# remotes::install_github("smeekes/bootUR")
# install.packages("forecast")
library(forecast)
library(bootUR)
library(tidyverse)
library(dplyr)
library(forecast)
library(ggplot2)

load("citibike.RData")

# --- Data Extraction -------

# Create a list
bike_ts_list <- list()

# 4. Loop through areas 1 to 6 to convert each into a 'ts' object
# Frequency = 24 captures the daily hourly pattern
for (i in 1:6) {
  col_name <- paste0("area ", i)

  bike_ts_list[[i]] <- ts(
    data_citibike[[col_name]],
    start = c(1, 1), # Starts at Day 1, Hour 1, which is hour 0 in our data
    frequency = 24 # 24 observations per day
  )
}

# 5. Extract individual time series variables for ease of use
ts_area1 <- bike_ts_list[[1]]
ts_area2 <- bike_ts_list[[2]]
ts_area3 <- bike_ts_list[[3]]
ts_area4 <- bike_ts_list[[4]]
ts_area5 <- bike_ts_list[[5]]
ts_area6 <- bike_ts_list[[6]]

# --- Data Exploration --------


for (i in 1:6) {
  png(paste0("outputs/hourly_demand_area_", i, ".png"),
    width = 1200, height = 800, res = 150
  )
  plot(bike_ts_list[[i]],
    main = paste("Hourly Bike Demand - Area", i),
    ylab = "Demand", xlab = "Time (Days)"
  )
  dev.off()
}


# visible: daily patterns
# the realized demand is only the people that actually got to rent the bike.
# The people that wanted to (have a demand) but for some reason couldnt, are not
# included in these numbers
# The true demand is 'latent', it is not observable.

# --- Unit root test ---------

packageVersion("bootUR")

area6 <- as.matrix(as.numeric(bike_ts_list[[6]]))

dim(area6)
class(area6)

d1_area6 <- diff(area6)
d2_area6 <- diff(d1_area6)

png("outputs/Hourly_demand_diffed.png", width = 1200, height = 800, res = 150)
plot(d1_area6,
  type = "l",
  main = "Hourly bike demand - once differentiated",
  xlab = "Time", ylab = "Demand",
  col = "magenta", lwd = 1.5
)
abline(h = 0, lty = 2, col = "lightblue")
dev.off()

png("outputs/Hourly_demand_twice_diffed.png", width = 1200, height = 800, res = 150)
plot(d2_area6,
  type = "l",
  main = "Hourly bike demand - twice differentiated",
  xlab = "Time", ylab = "Demand",
  col = "magenta", lwd = 1.5
)
abline(h = 0, lty = 2, col = "lightblue")
dev.off()


# Build predictors for the full sample --
df_full <- data_citibike
df_full$datetime <- ymd_h(paste(df_full$year, df_full$month, df_full$day, df_full$hour))
df_full$weekday <- wday(df_full$datetime, label = TRUE)
df_full$is_weekend <- factor(ifelse(df_full$weekday %in% c("Sat", "Sun"), "Weekend", "Workday"))

# use ALL levels from the full sample so train/test share the same factor levels
df_full$hour_f <- factor(df_full$hour)
df_full$month_f <- relevel(factor(df_full$month), ref = "4")

# attach area6 as a column so it splits together with the rest
df_full$area6 <- as.numeric(bike_ts_list[[6]])

# Same train/test split as before --
train_full <- df_full %>% dplyr::filter(datetime < ymd("2023-08-01"))
test_full <- df_full %>% dplyr::filter(datetime >= ymd("2023-08-01"))


# First differenced
train_full <- train_full[order(train_full$datetime), ]
area6_d1_train <- diff(train_full$area6)
is_weekend_d1 <- train_full$is_weekend[-1]
hour_f_d1 <- train_full$hour_f[-1]
month_f_d1 <- train_full$month_f[-1]

fit_d1_seasonal <- lm(area6_d1_train ~ 0 + is_weekend_d1:hour_f_d1 + month_f_d1)
d1_residuals <- residuals(fit_d1_seasonal)
test_d1 <- boot_adf(d1_residuals, deterministics = "intercept")
print(test_d1)

# Non differenced
fit_d0_seasonal <- lm(area6 ~ 0 + is_weekend:hour_f + month_f, data = train_full)
d0_residuals <- residuals(fit_d0_seasonal)
test_d0 <- boot_adf(d0_residuals, deterministics = "intercept")
print(test_d0)

# --- Exponential Smoothing Model ----------

y_train <- ts(train_full$"area 6", frequency = 24)
fit_ets <- ets(y_train, model = "ZZZ")
summary(fit_ets)

fit_ets$method
fit_ets$par


png("outputs/fit_ets", width = 1200, height = 800, res = 150)
plot(fit_ets)
dev.off()

# Ljung Box test
LB_test <- Box.test(residuals(fit_ets), lag = 24, type = "Ljung-Box")
LB_test

# --- Rolling Window ESM --------


# ---- Setup
y_all <- c(train_full$"area 6", test_full$"area 6") # full series
n_train <- length(train_full$"area 6")
H <- 24
n_days <- floor((length(y_all) - n_train) / H) # number of rolling origins

# Storage: one row per origin, one column per horizon
err_ets <- matrix(NA, n_days, H)
ets_forecast <- matrix(NA, n_days, H)

# ---- Rolling window -
for (d in 0:(n_days - 1)) {
  end <- n_train + d * H # last observed index at this origin
  y_win <- ts(y_all[(end - n_train + 1):end], frequency = 24) # rolling window
  actual <- y_all[(end + 1):(end + H)]

  # ETS: reuse the selected model form, re-estimate on the new window
  fit_d <- ets(y_win, model = fit_ets)
  f_ets <- forecast(fit_d, h = H)$mean

  err_ets[d + 1, ] <- actual - as.numeric(f_ets)
  ets_forecast[d + 1, ] <- as.numeric(f_ets)
}

# PUt all data points in 24h means
test_area6 <- test_full$area6
actual_matrix <- matrix(test_area6, nrow = 31, ncol = 24, byrow = TRUE)
actual_per_horizon <- colMeans(actual_matrix)
ets_f_per_horizon <- colMeans(ets_forecast)

results <- data.frame(
  h = 1:24,
  ETS = ets_f_per_horizon,
  actual = actual_per_horizon
)
print(results)

png("outputs/Demand_comp.png", width = 1200, height = 800, res = 150)
matplot(results$h, results[, c("ETS", "actual")],
  type = "l", lty = 1,
  col = c("magenta", "turquoise"), xlab = "Horizon h", ylab = "Demand"
)
legend("topleft", c("ETS Forecast", "Actual"), col = c("magenta", "turquoise"), lty = 1)
dev.off()


ets_for_vec <- as.vector(t(ets_forecast))

df_plot <- data.frame(
  datetime = seq(as.POSIXct("2023-08-01 00:00:00", tz = "UTC"),
    by = "hour", length.out = length(ets_for_vec)
  ),
  Actual = as.numeric(test_area6),
  Forecast = ets_for_vec
)

plt <- ggplot(df_plot, aes(x = datetime)) +
  geom_line(aes(y = Actual, colour = "Actual")) +
  geom_line(aes(y = Forecast, colour = "ETS rolling forecast")) +
  labs(x = NULL, y = "Demand", colour = NULL) +
  theme_minimal()

ggsave(
  plt,
  filename = file.path("outputs", "ets_vs_actual.png"),
  width = 10, height = 6
)


# Save results to further use when creating SARIMA
MAE_ETS <- (abs(err_ets))
MAE_p_horizon <- colMeans(MAE_ETS)
save(MAE_p_horizon, file = "outputs/MAE_p_horizon_ETS.RData")
err_ETS <- as.vector(t(err_ets))
save(err_ETS, file = "outputs/err_ETS.RData")
