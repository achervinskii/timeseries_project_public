library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
options(lubridate.week.start = 1)
load("./citibike.RData")
df <- data_citibike
head(df)

# pivot areas
df_long <- df |>
  mutate(datetime = make_datetime(year, month, day, hour)) |>
  pivot_longer(
    cols = paste0("area ", 1:6),
    names_to = "area",
    values_to = "demand",
    names_prefix = "area ",
  ) |>
  mutate(area = as.factor(area))
df_long

# compute aggregates on different horizons, keeping the area dimension
# mean hourly demand per month
df_month <- df_long |>
  mutate(
    month = month(datetime, label = TRUE),
  ) |>
  group_by(area, month) |>
  summarise(
    mean_demand = mean(demand),
    .groups = "drop"
  )

# mean daily demand per weekday
df_weekday <- df_long |>
  mutate(
    wday = wday(datetime,
      label = TRUE
    ),
  ) |>
  group_by(area, wday) |>
  summarise(mean_demand = mean(demand), .groups = "drop")


df_dayhour <- df_long |>
  mutate(hour = hour(datetime)) |>
  group_by(area, hour) |>
  summarise(
    mean_demand = mean(demand),
    .groups = "drop"
  )

plot_month <- ggplot(
  df_month,
  aes(x = month, y = mean_demand, group = area)
) +
  geom_line() +
  facet_wrap(~area, scales = "free_y") +
  labs(
    title = "Mean demand by month",
    x = "Month",
    y = "Mean hourly demand"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

plot_weekday <- ggplot(
  df_weekday,
  aes(x = wday, y = mean_demand, group = area)
) +
  geom_line() +
  facet_wrap(~area, scales = "free_y") +
  labs(
    title = "Mean demand by day of week",
    x = "Day of week",
    y = "Mean hourly demand"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

plot_hour <- ggplot(
  df_dayhour,
  aes(x = hour, y = mean_demand, group = area)
) +
  geom_line() +
  facet_wrap(~area, scales = "free_y") +
  labs(
    title = "Mean demand by hour of day",
    x = "Hour of day",
    y = "Mean hourly demand"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

# there is a pattern in the hour of day demand across different areas
# there is a pattern in the weekday, but it is lest pronounced across areas
# general upward trend over the months across different areas
# may be good to use ARIMA with dummies (weekdays and dayhours).

# look for interaction of a weekday and hour of the day
df_dayhour_weekday <- df_long |>
  filter(area == 6) |>
  mutate(weekday = wday(datetime, label = TRUE), hour = hour(datetime)) |>
  group_by(area, weekday, hour) |>
  summarise(mean_demand = mean(demand), .groups = "drop")

plot_hour_dayweek <- ggplot(
  df_dayhour_weekday,
  aes(weekday, hour, fill = mean_demand)
) +
  geom_tile() +
  labs(
    title = "Mean demand by day of week and hour of day, for area 6",
    x = "Day of week",
    y = "Hour of day"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

# the plot suggests a different pattern for weekends than for workdays
# save plots
plots <- list(
  plot_month = plot_month,
  plot_wday = plot_weekday,
  plot_hour = plot_hour,
  plot_hour_dayweek = plot_hour_dayweek
)

for (name in names(plots)) {
  ggsave(
    plots[[name]],
    filename = file.path("outputs", paste(name, ".png", sep = "")),
    width = 10, height = 6
  )
}


for (i in 1:6) {
  area_data <- df_long |>
    filter(area == i) |>
    select(demand, datetime)
  plt <- ggplot(area_data, aes(datetime, demand)) +
    geom_line() +
    labs(
      title = paste("Hourly Bike Demand - Area", i),
      y = "Demand", x = "Time (Days)"
    )
  name <- paste("plot_hourly_area_", i, sep = "")
  ggsave(
    plt,
    filename = file.path("outputs", paste(name, ".png", sep = "")),
    width = 10, height = 6
  )
}
# visible: daily patterns
# the realized demand is only the people that actually got to rent the bike.
# The people that wanted to (have a demand) but for some reason couldnt, are not
# included in these numbers
# The true demand is 'latent', it is not observable.
