library(sf)
library(sfdep)
library(spmodel)
library(ggplot2)
library(dplyr)

# Read the shapefile
dat <- st_read("HoustonHeat/Data/HoustonHeat.shp")

# Fix geometry issues (as noted in assignment)
dat <- st_make_valid(dat)


# Calculate overall rate R
R <- sum(dat$Count) / sum(dat$Population)

dat <- dat %>%
  mutate(
    EV     = Population * R,
    RR     = (Count + 0.5) / (EV + 0.5),
    log_RR = log(RR)
  )

# Density plot of log(RR)
ggplot(dat, aes(x = log_RR)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  labs(title = "Density of log(Relative Risk)",
       x = "log(RR)", y = "Density") +
  theme_minimal()



ggplot(dat, aes(fill = log_RR)) +
  geom_sf(color = NA) +
  scale_fill_viridis_c(option = "magma", name = "log(RR)") +
  labs(title = "Observed log(Relative Risk) of Heat-Related Illness",
       subtitle = "Houston Census Block Groups") +
  theme_void()


# Fit a basic OLS model
lm_fit <- lm(log_RR ~ NOAC + MED_AGE + HispanicPC + BlackPCT +
               under5PCT + over65PCT + povertyPCT + alonePCT + MinTemp,
             data = dat)

dat$lm_resid <- residuals(lm_fit)

# Chloropleth map of residuals
ggplot(dat, aes(fill = lm_resid)) +
  geom_sf(color = NA) +
  scale_fill_viridis_c(option = "plasma", name = "Residual") +
  labs(title = "OLS Residuals — Checking for Spatial Pattern") +
  theme_void()

# Moran's I test on residuals
nb <- st_contiguity(dat)
wt <- st_weights(nb)

global_moran_test(x  = dat$lm_resid,
                  nb = nb,
                  wt = wt)



car_fit <- spautor(
  log_RR ~ NOAC + MED_AGE + HispanicPC + BlackPCT +
    under5PCT + over65PCT + povertyPCT + alonePCT + MinTemp,
  data      = dat,
  spcov_type = "car"
)

summary(car_fit)
coef(car_fit)



# Standardized residuals
dat$std_resid <- rstandard(car_fit)

# 1. Linearity + Independence: Residuals vs Fitted
plot(fitted(car_fit), dat$std_resid,
     xlab = "Fitted Values", ylab = "Standardized Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)

# 2. Normality: Q-Q plot
qqnorm(dat$std_resid, main = "Normal Q-Q Plot")
qqline(dat$std_resid, col = "red")

# 3. Equal variance: similar spread across fitted values (check plot above)

# Model fit
pseudoR2(car_fit)


dat$fitted_RR <- fitted(car_fit)

ggplot(dat, aes(fill = fitted_RR)) +
  geom_sf(color = NA) +
  scale_fill_viridis_c(option = "magma", name = "Fitted log(RR)") +
  labs(title = "Model-Based Heat Illness Risk Across Houston",
       subtitle = "CAR Model Predictions by Census Block Group",
       caption = "Higher values = higher risk relative to expected") +
  theme_void()

# Save the fitted risk map as an image
ggsave("houston_risk_map.png",
       ggplot(dat, aes(fill = fitted_RR)) +
         geom_sf(color = NA) +
         scale_fill_viridis_c(option = "magma", name = "Fitted log(RR)") +
         labs(title = "Heat Illness Risk Across Houston") +
         theme_void(),
       width = 8, height = 6, dpi = 300
)
