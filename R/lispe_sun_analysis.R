## Code for White & Latty 2020

# Reset
  rm(list = ls())
  
# Libraries
  #install.packages(c('pavo', 'tidyverse', 'nlme', 'car', 'piecewiseSEM', 'mgcv', 'rstatix', 'ggpubr', 'inflection'))
  library(pavo)
  library(tidyverse)
  library(nlme)
  library(piecewiseSEM)
  library(car)
  library(mgcv)
  library(rstatix)
  library(ggpubr)
  library(inflection)
  
# Load data
  orient <- read.csv('../data/orientations.csv', stringsAsFactors = FALSE)  # fly orientations
  cloud <- procspec(as.rspec(read.csv('../data/spec_cloud.csv'), lim = c(300, 700)), opt = 'smooth', fixneg = 'zero')  # sunny irradiance data
  sun <- procspec(as.rspec(read.csv('../data/spec_sun.csv'), lim = c(300, 700)), opt = 'smooth', fixneg = 'zero')  # cloudy irradiance data
  lispe <- procspec(as.rspec(read.csv('../data/spec_lispe.csv'), lim = c(300, 700)), opt = 'smooth', fixneg = 'zero')  # aggregated fly reflectances
  
# Data processing
  lispe_face <- subset(lispe, 'face')  # subset fly faces
  lispe_wing <- subset(lispe, 'wing')  # subset fly wings
  lispe_face[, 2] <- lispe_face[, 2]/100  # normalise spectra
  lispe_wing[, 2] <- lispe_wing[, 2]/100  # normalise spectra
  orient$date_posix <- as.POSIXct(orient$date_posix)  # convert dates to...dates
  
  # Convert some data to long format for ease of plotting later
  sun_long <- gather(sun, orientation, reflectance, orient_000:orient_180, factor_key = TRUE)
  cloud_long <- gather(cloud, orientation, reflectance, orient_000:orient_180, factor_key = TRUE)
  lispe_long <- merge(lispe_face, lispe_wing)
  lispe_long <- gather(lispe_long, region, reflectance, face:wing, factor_key = TRUE)
  
# Subset data by treatment for convenient plotting
  orient_sun <- subset(orient, condition == 'sun')  # observational - full sun
  orient_cloud <- subset(orient, condition == 'cloud')  # observational - heavy cloud
  orient_sunblock <- subset(orient, condition == 'sunblock')  # manipulative - obscured sunlight
  
# Summarise orientation data
  orient_summ <- summarise(group_by(orient, condition),
                           n = n(),
                           n_days = length(unique(date)),
                           mean_offset = mean(offset),
                           median_offset = median(offset),
                           sd_offset = sd(offset),
                           se_offset = sd(offset)/sqrt(n),
                           min_solarexp = min(daily_solarexp),
                           max_solarexp = max(daily_solarexp),
                           mean_solarexp = mean(daily_solarexp),
                           sd_solarexp = sd(daily_solarexp),
                           min_rainfall = min(daily_rainfall),
                           max_rainfall = max(daily_rainfall),
                           mean_rainfall = mean(daily_rainfall),
                           sd_rainfall = sd(daily_rainfall),
                           min_maxtemp = min(daily_maxtemp),
                           max_maxtemp = max(daily_maxtemp),
                           mean_maxtemp = mean(daily_maxtemp),
                           sd_maxtemp = sd(daily_maxtemp))
  orient_summ
    
# Initial visualisation
  
  # Reorder orientation data relative to 'sun' condition for ease of plotting
  orient$condition <- factor(orient$condition, levels = c('sun', 'cloud', 'sunblock'))
  
  # Lispe face & wing reflectance
  ggplot(lispe_long, aes(wl, reflectance, group = region, colour = region)) +
    geom_line() +
    scale_x_continuous(breaks = seq(300, 700, 100)) +
    scale_y_continuous(breaks = seq(0, 1, 0.25)) +
    ylab('Reflectance (%)') +
    xlab('Wavelength (nm)') +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          panel.grid.major = element_line(),
          panel.grid.minor = element_blank(),
          legend.position = "none")
  
  ggsave('../figs/lispe_reflectance.tiff')
  ggsave('../figs/lispe_reflectance.png')
  
  # Male display offset ~ experimental condition
  ggplot(orient, aes(x = condition, y = offset)) +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    geom_jitter(width = 0.3, colour = 'darkgrey', alpha = 0.4) +
    scale_y_continuous(breaks = c(seq(0, 180, 45))) +
    ylab('Angular offset from solar azimuth (deg)') + 
    theme_bw() +
    theme(panel.grid.major = element_line(),
          panel.grid.minor = element_blank())
  
  ggsave('../figs/orient_condition.tiff')
  ggsave('../figs/orient_condition.png')
  
  # Male display offset ~ hour x experimental condition
  ggplot(orient, aes(x = hour, y = offset, group = hour)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(width = 15, colour = 'darkgrey', alpha = 0.4) +
    scale_y_continuous(breaks = seq(0, 180, 45)) +
    scale_x_continuous(breaks = seq(800, 1600, 100)) +
    ylab('Absolute offset from solar azimuth (deg)') + 
    #geom_hline(yintercept = 90, linetype = 'dashed') + 
    facet_grid(rows = vars(condition)) +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          panel.grid.major = element_line(),
          panel.grid.minor = element_blank())
  
  ggsave('../figs/orient_hour.tiff', height = 12)
  ggsave('../figs/orient_hour.png', height = 12)
 
# Estimate fly facial and wing luminance as a function of solar orientations across 180-degree range of orientations
  # Function for running the simple visual model
  colmod <- function(spec, sky){
    model <- lapply(2:20, function(x) vismodel(spec, visual = 'musca', achromatic = 'md.r1', illum = sky[, x]))
    space <- lapply(seq_along(model), function(x) colspace(model[[x]], space = 'categorical'))
    space
  }
  
  # Visually model wings & faces in cloud & sun
  lispe_face_sun <- colmod(lispe_face, sun)
  lispe_wing_sun <- colmod(lispe_wing, sun)
  lispe_face_cloud <- colmod(lispe_face, cloud)
  lispe_wing_cloud <- colmod(lispe_wing, cloud)

  # Rearrange the modelling data for ease of analysis/plotting
  lispe_sun <- data.frame(orient = seq(0, 180, 10),
                          do.call(rbind, lispe_face_sun)['lum'],
                          do.call(rbind, lispe_wing_sun)['lum'])
  names(lispe_sun) <- c('orient', 'face_lum', 'wing_lum')
  lispe_cloud <- data.frame(orient = seq(0, 180, 10),
                          do.call(rbind, lispe_face_cloud)['lum'],
                          do.call(rbind, lispe_wing_cloud)['lum'])
  names(lispe_cloud) <- c('orient', 'face_lum', 'wing_lum')
  lispe_sun_long <- gather(lispe_sun, region, lum, face_lum:wing_lum, factor_key = TRUE)
  lispe_cloud_long <- gather(lispe_cloud, region, lum, face_lum:wing_lum, factor_key = TRUE)
  
  # Run GAMs to model male signal luminance ~ solar offset 
  gam_face_sun <- gam(face_lum ~ s(orient), data = lispe_sun, method = 'REML')
  summary(gam_face_sun)
  gam_wing_sun <- gam(wing_lum ~ s(orient), data = lispe_sun, method = 'REML')
  summary(gam_wing_sun)
  gam_face_cloud <- gam(face_lum ~ s(orient), data = lispe_cloud, method = 'REML')
  summary(gam_face_cloud)  # n.s
  gam_wing_cloud <- gam(wing_lum ~ s(orient), data = lispe_cloud, method = 'REML')
  summary(gam_wing_cloud)  # n.s
  
  # Estimate inflection points for full-sun condition alone (since there are no inflection points under cloud)
  as.data.frame(edeci(lispe_sun$orient, gam_face_sun$fitted.values, index = 1, k = 5))
  as.data.frame(edeci(lispe_sun$orient, gam_wing_sun$fitted.values, index = 1, k = 5))
  
# Plots
  
  # Irradiance in sun
  a <- ggplot(sun_long, aes(wl, reflectance, group = orientation)) +
    geom_line() +
    scale_x_continuous(breaks = seq(300, 700, 100)) +
    ylab('Absolute irradiance (μmol/s-1/m-2)') +
    xlab('Wavelength (nm)') +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          panel.grid.major = element_line(),
          panel.grid.minor = element_blank(),
          legend.position = "none")
  
  # Irradiance in cloud
  b <- ggplot(cloud_long, aes(wl, reflectance, group = orientation)) +
    geom_line() +
    scale_x_continuous(breaks = seq(300, 700, 100)) +
    ylim(c(0, 0.3)) +
    ylab(' ') +
    xlab('Wavelength (nm)') +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          panel.grid.major = element_line(),
          panel.grid.minor = element_blank(),
          legend.position = "none")

  # Signal luminance in the sun
  c <- ggplot(lispe_sun_long, aes(orient, lum, group = region, colour = region)) +
    geom_point() +
    geom_smooth(method = gam, formula = y ~ s(x)) +
    scale_x_continuous(breaks = seq(0, 180, 20)) +
    ylab('Signal luminance') +
    xlab('Absolute offset from solar azimuth (deg)') +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          panel.grid.major = element_line(),
          panel.grid.minor = element_blank(),
          legend.position = "none")
  
  # Signal luminance in cloud
  d <- ggplot(lispe_cloud_long, aes(orient, lum, group = region, colour = region)) +
    geom_point() +
    geom_smooth(method = gam, formula = y ~ s(x)) +
    scale_x_continuous(breaks = seq(0, 180, 20)) +
    ylab(' ') +
    xlab('Absolute offset from solar azimuth (deg)') +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          panel.grid.major = element_line(),
          panel.grid.minor = element_blank(),
          legend.position = "none")
  
  # Arrange & save
  figure <- ggarrange(a, b, c, d,
                      labels = c("A", "B", "C", "D"),
                      ncol = 2, nrow = 2)
  ggsave('../figs/spectra.tiff', width = 8, plot = figure)
  ggsave('../figs/spectra.png', width = 8, plot = figure)

## Density plots for melding into the orientation figure (fig 1).
 # Can't be bothered working out how to do the whole thing in R.
  ggplot(orient_sun, aes(offset)) +
    geom_histogram(color="darkblue", fill="lightblue") +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position = "none")
  ggsave('../figs/hist_sun.tiff', width = 15)
  
  ggplot(orient_cloud, aes(offset)) +
    geom_histogram(color="darkblue", fill="lightblue") +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position = "none")
  ggsave('../figs/hist_cloud.tiff', width = 15)

## Statistics
  
  # Use sunny GAMs to predict realised signal luminance in across actual male orientations
  orient$face_lum_sun <- as.numeric(predict(gam_face_sun, data.frame(orient = orient$offset)))
  orient$wing_lum_sun <- as.numeric(predict(gam_wing_sun, data.frame(orient = orient$offset)))
  
  # Reorder relative to cloud for ease of testing
  orient$condition <- factor(orient$condition, levels = c('cloud', 'sun', 'sunblock'))
  
  ## Main model of orientation ~ treatment. Can't use gamma dist because there are real 0's in the data, 
  # nor can use Gaussian with link = log for the same reason, so Gaussian with sqrt transform it is.
  m1 <- lme(sqrt(offset) ~ condition * sun_al, random = ~1 | date, data = orient)
  
  # Inspect the residuals. Looks alright.
  qqnorm(resid(m1))
  plot(m1)
  
  # Likelihood-ratio Chi^2 on overall model
  car::Anova(m1)
  
  # Model summary
  summary(m1)
  rsquared(m1)
  
  # Save output
  sink("../output/m1.txt")
  car::Anova(m1)
  summary(m1)
  rsquared(m1)
  sink()
  
  ## Test signal luminance ~ behaviour in sun only, given the lack of relationship between orientation & luminance under cloud
  orient_sub <- subset(orient, condition != 'sunblock')
  kruskal.test(face_lum_sun ~ condition, data = orient_sub)
  kruskal.test(wing_lum_sun ~ condition, data = orient_sub)
  
  ## Effect sizes
  # ETA squared
  kruskal_effsize(face_lum_sun ~ condition, data = orient_sub, ci = TRUE)
  kruskal_effsize(wing_lum_sun ~ condition, data = orient_sub, ci = TRUE)
  
  # Luminance difference
  median(subset(orient_sub, condition == 'sun')[['face_lum_sun']]) / median(subset(orient_sub, condition == 'cloud')[['face_lum_sun']])
  median(subset(orient_sub, condition == 'sun')[['wing_lum_sun']]) / median(subset(orient_sub, condition == 'cloud')[['wing_lum_sun']])

  