
library(mgcv); library(data.table); library(gratia); library(marginaleffects); library(ggplot2); library(dplyr); library(tidyr); library(zoo); library(slider)



#####test Paracou####

obs <- read.csv("G:/My Drive/cesab bioforest/Data/aggregated_data_v9.csv")

setDT(obs)

obs=obs[obs$Site=="Paracou"]

unique(obs$Plot)

setnames(obs, "Plot", "subplot")

obs$plot=gsub("_.", "",obs$subplot )

obs <- obs[!Year%in%1984:1990]

obs <- dcast(
  obs,
  subplot + plot + Year ~ variable,
  value.var = "value"
)


##### upload data 
load("clim_by_site.rda")

clim_all=rbindlist(lapply(clim_by_site[["Paracou"]], function(x) x$climate[[1]]),
                   idcol = "plot")
#which period to consider the anomalies? Here I take the from 1977 (the oldest inventory date among all sites considered)

#"Paracou", "Mbaiki", "Corinto", "SUAS", "Lesong", "Ulu Muda", "Sungai Lalang", "Tene", "Jari"

clim_stats <- clim_all[, .(
  mean_tmax = mean(tmax),
  sd_tmax   = sd(tmax),
  mean_vpd   = mean(vpd),
  sd_vpd     = sd(vpd),
  mean_srad  = mean(srad),
  sd_srad    = sd(srad),
  mean_def  = mean(def),
  sd_def    = sd(def)
), by = .(plot, month)]

clim_all <- merge(clim_all, clim_stats, by = c("plot", "month"), all.x = TRUE)

clim_all[, `:=`(
  z_tmax = (tmax - mean_tmax) / sd_tmax,
  z_vpd  = (vpd  - mean_vpd)  / sd_vpd,
  z_srad = (srad - mean_srad) / sd_srad,
  z_def  = (def  - mean_def)  / sd_def
)]

#climate anomalies follow inventories. For Paracou, mean per year before 1996, two years after

clim_all[, year_bin :=
           ifelse(year <= 1995,
                  year,                     # keep original year
                  1997 + ((year - 1996) %/% 2) * 2   # 2-year bins
           )]


census_anom <- clim_all[,.(
  mean_z_tmax = mean(z_tmax, na.rm = TRUE),
  mean_z_vpd = mean(z_vpd, na.rm = TRUE),
  mean_z_srad = mean(z_srad, na.rm = TRUE),
  mean_z_def = mean(z_def, na.rm = TRUE)
),
by = .(plot, year_bin)]

#AGB values before 1991 are not consistent
census_anom=census_anom[year_bin>=1991] 

unique(census_anom$plot)
unique(obs$plot)

merged <- merge(
  obs,
  census_anom,
  by.x = c("Year","plot"), by.y=c("year_bin","plot")
)


control=c(1,6,11,13,14,15)



merged$Treatment <- ifelse(merged$plot %in% control, "control", "disturbed")

control=merged[Treatment=="control",]

disturbed=merged[Treatment=="disturbed",]

merged$Treatment <- factor(merged$Treatment)
disturbed$plot <- factor(disturbed$plot)
disturbed$subplot <- factor(disturbed$subplot)


######fit trajectory only
####GAM



disturbed <- disturbed %>%
  group_by(subplot) %>%
  mutate(
    fitted_GAM5 = predict(
      gam(agb ~ s(Year, k = 5), data = cur_data()),
      newdata = cur_data()
    )
  ) %>%
  ungroup()


disturbed <- disturbed %>%
  group_by(subplot) %>%
  mutate(
    fitted_GAM8 = predict(
      gam(agb ~ s(Year, k = 8), data = cur_data()),
      newdata = cur_data()
    )
  )%>%
  ungroup()



p1=ggplot(disturbed, aes(x = Year, y = fitted_GAM5, color = subplot)) +
  geom_line(size = 1) +  # predicted trajectory
  geom_point(data = disturbed, 
             aes(x = Year, y = agb, color = subplot),
             inherit.aes = FALSE, size = 2) +  # observed points
  theme_bw() +
  labs(x = "Year", 
       y = "AGB", title="GAM-k=5") + theme(legend.position = "none")


plot(disturbed$agb, disturbed$fitted_GAM5)

p1.1 <- ggplot(disturbed, aes(Year, agb - fitted_GAM5)) +
  geom_point() +
  ggtitle("GAM residuals")

p1 + p1.1


p11=ggplot(disturbed, aes(x = Year, y = fitted_GAM8, color = subplot)) +
  geom_line(size = 1) +  # predicted trajectory
  geom_point(data = disturbed, 
             aes(x = Year, y = agb, color = subplot),
             inherit.aes = FALSE, size = 2) +  # observed points
  theme_bw() +
  labs(x = "Year", 
       y = "AGB", title="GAM-k=8") + theme(legend.position = "none")


plot(disturbed$agb, disturbed$fitted_GAM8)

p11.1 <- ggplot(disturbed, aes(Year, agb - fitted_GAM8)) +
  geom_point() +
  ggtitle("GAM residuals")

p11 + p11.1



#####rolling average

#to remove NAs

disturbed <- disturbed %>%
  group_by(subplot) %>%
  arrange(Year) %>%
  mutate(fitted_rollmean2 = zoo::rollapply(
    agb,
    width = 2,
    FUN = mean,
    align = "center",
    fill = NA,
    partial = TRUE  #uses data to avoid NA
  ))%>%
  ungroup()

disturbed <- disturbed %>%
  group_by(subplot) %>%
  arrange(Year) %>%
  mutate(fitted_rollmean5 = zoo::rollapply(
    agb,
    width = 5,
    FUN = mean,
    align = "center",
    fill = NA,
    partial = TRUE  #
  ))%>%
  ungroup()



p2=ggplot(disturbed, aes(x = Year, y = fitted_rollmean2, color = subplot)) +
  geom_line(size = 1) +  # predicted trajectory
  geom_point(data = disturbed, 
             aes(x = Year, y = agb, color = subplot),
             inherit.aes = FALSE, size = 2) +  # observed points
  theme_bw() +
  labs(x = "Year", 
       y = "AGB", title="rollmean - k=2 align =centre") + theme(legend.position = "none")

plot(disturbed$agb, disturbed$fitted_rollmean2)

p2.1 <- ggplot(disturbed, aes(Year, agb - fitted_rollmean2)) +
  geom_point() +
  ggtitle("rollmean residuals")

p2 + p2.1


p22=ggplot(disturbed, aes(x = Year, y = fitted_rollmean5, color = subplot)) +
  geom_line(size = 1) +  # predicted trajectory
  geom_point(data = disturbed, 
             aes(x = Year, y = agb, color = subplot),
             inherit.aes = FALSE, size = 2) +  # observed points
  theme_bw() +
  labs(x = "Year", 
       y = "AGB", title="rollmean - k=5 align =centre") + theme(legend.position = "none")

plot(disturbed$agb, disturbed$fitted_rollmean5)

p22.1 <- ggplot(disturbed, aes(Year, agb - fitted_rollmean5)) +
  geom_point() +
  ggtitle("rollmean residuals")

p22 + p22.1


###slider


disturbed <- disturbed %>%
  group_by(subplot) %>%
  arrange(Year) %>%
  mutate(fitted_slide = slide_dbl(agb, mean, .before = 2, .after = 2))%>%
  ungroup()


p3=ggplot(disturbed, aes(x = Year, y = fitted_slide, color = subplot)) +
  geom_line(size = 1) +  # predicted trajectory
  geom_point(data = disturbed, 
             aes(x = Year, y = agb, color = subplot),
             inherit.aes = FALSE, size = 2) +  # observed points
  theme_bw() +
  labs(x = "Year", 
       y = "AGB", title="slide - .before = 2, .after = 2")+ theme(legend.position = "none")

plot(disturbed$agb, disturbed$fitted_slide)

p3.1 <- ggplot(disturbed, aes(Year, agb - fitted_slide)) +
  geom_point() +
  ggtitle("slide residuals")

p3 + p3.1


#polynome


disturbed <- disturbed %>%
  group_by(subplot) %>%
  mutate(fitted_poly2 = predict(lm(agb ~ poly(Year, 2), data = cur_data())))%>%
  ungroup()

p4=ggplot(disturbed, aes(x = Year, y = fitted_poly2, color = subplot)) +
  geom_line(size = 1) +  # predicted trajectory
  geom_point(data = disturbed, 
             aes(x = Year, y = agb, color = subplot),
             inherit.aes = FALSE, size = 2) +  # observed points
  theme_bw() +
  labs(x = "Year", 
       y = "AGB", title="Quad poly")+ theme(legend.position = "none")

plot(disturbed$agb, disturbed$fitted_poly2)

p4.1 <- ggplot(disturbed, aes(Year, agb - fitted_poly2)) +
  geom_point() +
  ggtitle("Poly residuals")

p4 + p4.1




#loess


disturbed <- disturbed %>%
  group_by(subplot) %>%
  mutate(fitted_LOESS0.75 = predict(loess(agb ~ Year,
                                   data = cur_data(),
                                   span = 0.75)))%>%
  ungroup()

disturbed <- disturbed %>%
  group_by(subplot) %>%
  mutate(fitted_LOESS0.5 = predict(loess(agb ~ Year,
                                          data = cur_data(),
                                          span = 0.5)))%>%
  ungroup()



p5=ggplot(disturbed, aes(x = Year, y = fitted_LOESS0.75, color = subplot)) +
  geom_line(size = 1) +  # predicted trajectory
  geom_point(data = disturbed, 
             aes(x = Year, y = agb, color = subplot),
             inherit.aes = FALSE, size = 2) +  # observed points
  theme_bw() +
  labs(x = "Year", 
       y = "AGB", title="LOESS - 0.75")+ theme(legend.position = "none")

plot(disturbed$agb, disturbed$fitted_LOESS0.75)

p5.1 <- ggplot(disturbed, aes(Year, agb - fitted_LOESS0.75)) +
  geom_point() +
  ggtitle("LOESS residuals")

p5 + p5.1

p55=ggplot(disturbed, aes(x = Year, y = fitted_LOESS0.5, color = subplot)) +
  geom_line(size = 1) +  # predicted trajectory
  geom_point(data = disturbed, 
             aes(x = Year, y = agb, color = subplot),
             inherit.aes = FALSE, size = 2) +  # observed points
  theme_bw() +
  labs(x = "Year", 
       y = "AGB", title="LOESS - 0.05")+ theme(legend.position = "none")

plot(disturbed$agb, disturbed$fitted_LOESS0.5)

p55.1 <- ggplot(disturbed, aes(Year, agb - fitted_LOESS0.5)) +
  geom_point() +
  ggtitle("LOESS residuals")

p55 + p55.1




names(disturbed)

#root mean square error

rmse_results <- disturbed %>%
  group_by(subplot) %>%
  summarise(
    RMSE_GAM8   = sqrt(mean((agb - fitted_GAM8)^2)),
    RMSE_GAM5   = sqrt(mean((agb - fitted_GAM5)^2)),
    RMSE_ROLLMEAN2 = sqrt(mean((agb - fitted_rollmean2)^2)),
    RMSE_ROLLMEAN5 = sqrt(mean((agb - fitted_rollmean5)^2)),
    RMSE_SLIDE  = sqrt(mean((agb - fitted_slide)^2)),
    RMSE_LOESS0.75  = sqrt(mean((agb - fitted_LOESS0.75)^2)),
    RMSE_LOESS0.5  = sqrt(mean((agb - fitted_LOESS0.5)^2)),
    RMSE_PolyQuad  = sqrt(mean((agb - fitted_poly2)^2))
  ) %>% 
  ungroup()

rmse_results <- rmse_results %>%
  mutate(
    best_method = colnames(.)[
      apply(select(., -subplot), 1, which.min) + 1
    ]
  )

as.data.table(rmse_results)

vars <- grep("^fitted_", names(disturbed), value = TRUE)

for (v in vars) {
  cat(v, "\n")
  print(summary(lm(as.formula(paste0(v, " - agb ~ Year")), data = disturbed)))
}


vars <- grep("^fitted_", names(disturbed), value = TRUE)

results <- data.frame(
  method = vars,
  intercept = NA,
  slope = NA,
  p_value = NA
)

for (i in seq_along(vars)) {
  
  v <- vars[i]
  model <- lm(as.formula(paste0(v, " - agb ~ Year")), data = disturbed)
  s <- summary(model)
  
  results$intercept[i] <- coef(s)[1, 1]
  results$slope[i]     <- coef(s)[2, 1]
  results$p_value[i]   <- coef(s)[2, 4]
}

results





pairwise_corr_for <- function(df, colname) {
  
  wide <- df %>%
    select(subplot, Year, value = all_of(colname)) %>%
    filter(!is.na(value)) %>%
    pivot_wider(names_from = subplot, values_from = value)
  
  M <- as.matrix(wide[ , -1])             # drop Year column
  
  C <- cor(M, use = "pairwise.complete.obs")
  
  mean(C[upper.tri(C)], na.rm = TRUE)     # average off-diagonal correlations
}


methods <- c(
  "fitted_GAM5",
  "fitted_GAM8",
  "fitted_rollmean2",
  "fitted_rollmean5",
  "fitted_slide",
  "fitted_LOESS0.75",
  "fitted_LOESS0.5",
  "fitted_poly2"
)

mean_corr <- sapply(methods, function(m) pairwise_corr_for(disturbed, m))
temporal_signal <- data.frame(
  method = methods,
  mean_pairwise_corr = mean_corr
)
temporal_signal <- temporal_signal[order(-temporal_signal$mean_pairwise_corr), ]

temporal_signal

#dispersion error dans le temps
#polynome
