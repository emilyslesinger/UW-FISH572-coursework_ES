##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       fit index model to the simulated data
## Author:        Primarily Lewis Barnett
## Notes:         Turn on 'outline' to see structure of code
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load packages ----------------------------------------------------------------
library(dplyr) # install.packages("dplyr")
library(ggplot2) # install.packages("ggplot2")
library(sdmTMB) # install.packages("sdmTMB")
library(here) # install.packages("here")

# Load simulated data ----------------------------------------------------------
sim_dat <- readRDS(here::here("coursework/simulations/sim_data/sim_dat_20.RDS"))

# sample 100 locations from each year using simple random sampling
set.seed(99)
sim_dat_obs <- sim_dat %>%
  dplyr::group_by(year) %>%
  dplyr::slice_sample(n = 100) 

# Explore simulated data -------------------------------------------------------
# plot samples over surface of mean without observation error (taken as "true")

ggplot2::ggplot(sim_dat, aes(X, Y)) +
  ggplot2::geom_raster(aes(fill = eta_scaled)) +
  ggplot2::geom_point(aes(size = observed_scaled), data = sim_dat_obs, pch = 21) +
  ggplot2::facet_wrap(~year) +
  ggplot2::scale_fill_viridis_c() +
  ggplot2::scale_size_area() +
  ggplot2::coord_cartesian(expand = FALSE)

# Fit new model ----------------------------------------------------------------
mesh <- sdmTMB::make_mesh(sim_dat_obs, 
                          xy_cols = c("X", "Y"), 
                          type = "cutoff_search", 
                          n_knots = 50)
plot(mesh)

fit <- sdmTMB::sdmTMB(
  formula = observed_scaled ~ 0 + as.factor(year), 
  data = sim_dat_obs,
  mesh = mesh,
  time = "year",
  family = sdmTMB::tweedie(), 
  spatial = "on", # c("on", "off")
  spatiotemporal = "iid", # c("iid", "ar1", "rw", "off")
)

fit

sanity(fit) # model checking

# Inspect model ----------------------------------------------------------------
# randomized quantile residuals
sim_dat_obs$resids <- residuals(fit) 
hist(sim_dat_obs$resids)

# qq plot
qqnorm(sim_dat_obs$resids)
abline(a = 0, b = 1)

# spatial residuals
ggplot(sim_dat_obs, aes(X, Y, col = resids)) + scale_colour_gradient2() +
  geom_point() + facet_wrap(~year) + coord_fixed()

# predict model to grid of full survey domain ----------------------------------
# replicate grid for each year to make prediction grid
grid <- readRDS(here::here("coursework/simulations/sim_data/grid.RDS"))
grid_yrs <- sdmTMB::replicate_df(grid, "year", unique(sim_dat$year))

# predict
predictions <- stats::predict(fit, newdata = grid_yrs, return_tmb_object = TRUE)

# visualize predictions, then how fixed and random effects contribute 
plot_map <- function(dat, column) {
  ggplot(dat, aes(X, Y, fill = {{ column }})) +
    geom_raster() +
    facet_wrap(~year) +
    coord_fixed()
}

# full prediction
plot_map(predictions$data, exp(est)) +
  scale_fill_viridis_c(trans = "sqrt") +
  ggtitle("Prediction (fixed effects + all random effects)")

# fixed effects only (year)
plot_map(predictions$data, exp(est_non_rf)) +
  ggtitle("Prediction (fixed effects only)") +
  scale_fill_viridis_c(trans = "sqrt")

# spatial random effects
plot_map(predictions$data, omega_s) +
  ggtitle("Spatial random effects only") +
  scale_fill_gradient2()

# spatiotemporal random effects
plot_map(predictions$data, epsilon_st) +
  ggtitle("Spatiotemporal random effects only") +
  scale_fill_gradient2()

# compute abundance index ------------------------------------------------------
# we will assume that the area of each grid cell is 1
index <- sdmTMB::get_index(predictions, area = 1, bias_correct = TRUE)

# plot index
ggplot2::ggplot(index, aes(year, est)) + 
  ggplot2::geom_line() +
  ggplot2::geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  ggplot2::xlab('Year') + 
  ggplot2::ylab('Total biomass estimate')
