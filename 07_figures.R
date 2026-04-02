#07_main figures and tables  

#load packages 
library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(RColorBrewer)
library(scales)

#load topo csv for gauge coordinates 
topo <- read.csv("Data/CAMELS/camels-20250724T0232Z/camels_topo.txt", sep = ";")
topo2 <- topo %>%
  mutate(id = sprintf("%08d", gauge_id))

#Background map 
#load usa and states (from naturalearth)
usa    <- ne_countries(scale = "medium", returnclass = "sf") |> subset(admin == "United States of America")
states <- ne_states(country = "united states of america", returnclass = "sf")

#Bounding Box (CONUS):
xlim <- c(-125, -66)
ylim <- c(24, 50)

##PLOTS

#01 CAMELS topograhic and climate attributes

#Elevation and slope 
#converts into sf_object; wgs84
topo_sf <- st_as_sf(topo, coords = c("gauge_lon", "gauge_lat"), crs = 4326)

#define colours
cols <- c(
  "#f5f5f5",
  "#f6e8c3",
  "#dfc27d",
  "#bf812d",
  "#8c510a",
  "#543005",
  "#1a0d00")

ggplot() +
  geom_sf(data = usa, fill = "grey95", color = "grey50", linewidth = 0.2) +
  geom_sf(data = states, fill = NA, color = "grey50", linewidth = 0.15) +
  geom_sf(
    data = topo_sf,
    aes(fill = elev_mean),
    shape = 21,
    colour = "grey35",
    stroke = 0.15,
    size = 2.8
  ) +
  scale_fill_stepsn(
    colours = cols,
    values  = rescale(c(0, 150, 300, 500, 1000, 2000, 3000, 4000), from = c(0, 4000)),
    breaks  = c(0, 150, 300, 500, 1000, 2000, 3000, 4000),
    limits  = c(0, 4000),
    labels  = c("0", "150", "300", "500", "1000", "2000", "3000", "4000"),
    name    = "Mean catchment elevation [m a.s.l.]",
    guide   = guide_coloursteps(
      title.position = "top",
      label.position = "bottom",
      even.steps = TRUE,
      barwidth = unit(8, "cm"),
      barheight = unit(0.5, "cm")
    )
  ) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.justification = "left",
    legend.box.just = "left",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

ggsave("Figures/CAMELS_elev.png", width = 8, height = 5, dpi = 600)

ggplot() +
  geom_sf(data = usa, fill = "grey95", color = "grey50", linewidth = 0.2) +
  geom_sf(data = states, fill = NA, color = "grey50", linewidth = 0.15) +
  geom_sf(
    data = topo_sf,
    aes(fill = slope_mean),
    shape = 21,
    colour = "grey35",
    stroke = 0.15,
    size = 2.8
  ) +
  scale_fill_stepsn(
    colours = cols,
    values  = rescale(c(0, 10, 25, 50, 100, 150, 200, 300), from = c(0, 300)),
    breaks  = c(0, 10, 25, 50, 100, 150, 200, 300),
    limits  = c(0, 300),
    labels  = c("0", "10", "25", "50", "100", "150", "200", "300"),
    name    = "Mean catchment slope [m km-1]",
    guide   = guide_coloursteps(
      title.position = "top",
      label.position = "bottom",
      even.steps = TRUE,
      barwidth = unit(8, "cm"),
      barheight = unit(0.5, "cm")
    )
  ) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.justification = "left",
    legend.box.just = "left",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

ggsave("Figures/CAMELS_slope.png", width = 8, height = 5, dpi = 600)

#Hydroclimatic classes

clim <- read.csv("Data/CAMELS/camels-20250724T0232Z/camels_clim.txt", sep = ";")

clim <- clim %>%
  mutate(id = sprintf("%08d", gauge_id))

#join with topo2 for coordinates 
clim_df <- clim %>%
  left_join(topo2, by = "id")

clim_sf <- st_as_sf(clim_df, coords = c("gauge_lon", "gauge_lat"), crs = 4326)

#define climate classes 
clim_cat <- clim_sf %>%
  mutate(
    hydro_class = case_when(
      frac_snow > 0.2 ~ "snow",
      aridity >= 1    ~ "arid",
      aridity < 1     ~ "humid"
    )
  )

ggplot() +
  geom_sf(data = usa, fill = "grey95", color = "grey50", linewidth = 0.2) +
  geom_sf(data = states, fill = NA, color = "grey50", linewidth = 0.15) +
  geom_sf(
    data = clim_cat,
    aes(fill = hydro_class),
    shape = 21,
    colour = "grey35",
    stroke = 0.15,
    size = 2.8, 
    alpha = 0.7
  ) +
  scale_fill_manual(
    values = c(
      snow  = "blue",
      arid  = "#F0E442",
      humid = "#66A61E"
    ),
    breaks = c("snow", "arid", "humid"),
    labels = c("snow", "arid", "humid"),
    name = "Hydroclimatic class"
  ) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.justification = "left",
    legend.box.just = "left",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

ggsave("Figures/hydro_class.png", width = 8, height = 5, dpi = 600)

#02 HBV model performance with CAMELS-forcing; CAMELS Changepoint analysis  
#HBV NSE (Camels)
nse <- read.csv("Output/CAMELS/Camels_cal/best_params_all_EZG_CAMELS.csv", header = TRUE)

nse <- nse %>%
  mutate(id = sprintf("%08d", id))

#join with topo2 for coordinates 
nse_points_df <- nse %>%
  left_join(topo2, by = "id")

nse_points_sf <- st_as_sf(nse_points_df, coords = c("gauge_lon", "gauge_lat"), crs = 4326)

ggplot() +
  geom_sf(data = usa, fill = "grey95", color = "grey50", linewidth = 0.2) +
  geom_sf(data = states, fill = NA, color = "grey50", linewidth = 0.15) +
  geom_sf(
    data = nse_points_sf,
    aes(fill = NSE),
    shape = 21,
    colour = "grey35",
    stroke = 0.15,
    size = 2.8
  ) +
  scale_fill_stepsn(
    colours = brewer.pal(9, "YlGnBu"),
    limits  = c(0, 1),
    breaks  = seq(0, 1, 0.2),
    labels  = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"),
    name    = "NSE (HBV with CAMELS)",
    guide   = guide_coloursteps(
      title.position = "top",
      label.position = "bottom",
      even.steps = TRUE,
      barwidth = unit(8, "cm"),
      barheight = unit(0.5, "cm")
    )
  ) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.justification = "left",
    legend.box.just = "left",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

ggsave("Figures/HBV_nse.png", width = 8, height = 5, dpi = 600)

#CAMELS Changepoint
CP_sum <- read.csv("Output/CAMELS/CAMELS_output/CP_CAMELS_full.csv")

#Changepoints
CP_sum <- CP_sum %>%
  mutate(id = sprintf("%08d", id))

CP_points_df <- CP_sum %>%
  left_join(topo2, by = "id")

CP_points_sf <- st_as_sf(CP_points_df, coords = c("gauge_lon", "gauge_lat"), crs = 4326)

my_blues <- c("#f7fbff", "#c6dbef","#6baed6", "#2171b5", "#08306b")

ggplot() +
  geom_sf(data = usa, fill = "grey95", color = "grey50", linewidth = 0.2) +
  geom_sf(data = states, fill = NA, color = "grey50", linewidth = 0.15) +
  geom_sf(data = CP_points_sf, aes(fill = CP), shape = 21, color = "grey35", size = 2.8, stroke = 0.15) +
  scale_fill_stepsn(
    colours = my_blues,
    limits  = c(0, 1),
    breaks  = seq(0, 1, 0.2),
    name    = "Changepoint (CAMELS)", 
    na.value = "grey70",
    guide   = guide_coloursteps(
      title.position = "top",
      label.position = "bottom",
      even.steps = TRUE,
      barwidth = unit(8, "cm"),
      barheight = unit(0.5, "cm"))
  ) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.justification = "left",
    legend.box.just = "left",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

ggsave("Figures/cp.png", width = 8, height = 5, dpi = 600)

#Changepoint nRMSE
ggplot() +
  geom_sf(data = usa, fill = "grey95", color = "grey50", linewidth = 0.2) +
  geom_sf(data = states, fill = NA, color = "grey50", linewidth = 0.15) +
  geom_sf(data = CP_points_sf, aes(fill = nRMSE), shape = 21, color = "grey35", size = 2.8, stroke = 0.15) +
  scale_fill_stepsn(
    colours = rev(brewer.pal(5, "YlGnBu")),
    limits  = c(0, 1),
    breaks  = seq(0, 1, 0.2),
    name    = "nRMSE Changepoint (CAMELS)", 
    na.value = "grey70",
    guide   = guide_coloursteps(
      title.position = "top",
      label.position = "bottom",
      even.steps = TRUE,
      barwidth = unit(8, "cm"),
      barheight = unit(0.5, "cm"))
  ) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.justification = "left",
    legend.box.just = "left",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

ggsave("Figures/cp_nrmse.png", width = 8, height = 5, dpi = 600)

#histogram changepoint distribution 
cp_med <- median(CP_sum$CP, na.rm = TRUE)

ggplot(CP_sum, aes(x = CP)) +
  geom_histogram(
    binwidth = 0.05,
    boundary = 0,
    fill = "grey40",
    color = "white"
  ) +
  geom_vline(xintercept = cp_med, linetype = "dashed", linewidth = 0.8) +
  annotate(
    "text",
    x = cp_med, y = Inf,
    label = paste0("Median = ", round(cp_med, 2)),
    vjust = 1.4, hjust = -0.1
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.1),
    minor_breaks = seq(0, 1, by = 0.05),   
    labels = number_format(accuracy = 0.1)
  ) +
  scale_y_continuous(
    limits = c(0, 125),
    breaks = seq(0, 120, by = 20),
    minor_breaks = seq(0, 120, by = 10),   
    expand = c(0, 0)
  ) +
  labs(
    x = "Soil moisture changepoint (-)",
    y = "Number of catchments"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major = element_line(color = "grey80", linewidth = 0.4)
  )
ggsave("Figures/CP_hist.png", width = 8, height = 6, dpi = 300)

##ISIMIP
#03 HBV model performance with ISIMIP-forcing 
nse <- read.csv("Output/ISIMIP/ISIMIP_cal/all_EZG_best_param_ISIMIP.csv", header = TRUE)

nse <- nse %>%
  mutate(id = sprintf("%08d", id))

#join with topo2 for coordinates 
nse_points_df <- nse %>%
  left_join(topo2, by = "id")

nse_points_sf <- st_as_sf(nse_points_df, coords = c("gauge_lon", "gauge_lat"), crs = 4326)

#categorize into nse categories
nse_points_sf <- nse_points_sf %>% 
  mutate(
    NSE_cat = cut(
      NSE,
      breaks = c(-Inf, -0.3, 0, 0.3, 0.6, 1),   
      labels = c("< -0,3", "-0.3-0" ,"0–0.3", "0.3–0.6", "> 0.6")
    )
  )

ggplot() +
  geom_sf(data = usa, fill = "grey95", color = "grey50", linewidth = 0.2) +
  geom_sf(data = states, fill = NA, color = "grey50", linewidth = 0.15) +
  geom_sf(
    data = nse_points_sf,
    aes(fill = NSE_cat),
    shape = 21, color = "grey35", size = 2.8, stroke = 0.15
  ) +
  scale_fill_brewer(
    palette = "YlGnBu",
    name    = "NSE (HBV with ISIMIP)"
  ) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.justification = "left",
    legend.box.just = "left",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

ggsave("Figures/HBV_nse_ISIMIP.png", width = 8, height = 5, dpi = 600)

#04 ISIMIP Data evaluation 
#04.1 ISIMIP seasonality : circular stats - mean occurence date 
cir_stat_sum <- read.csv("Output/ISIMIP/ISIMIP_output/cir_stat_sum.csv")

cir_stat_sum <- cir_stat_sum %>%
  mutate(id = sprintf("%08d", id))

cirstat_points_df <- cir_stat_sum %>%
  left_join(topo2, by = "id")

cirstat_points_sf <- st_as_sf(cirstat_points_df, coords = c("gauge_lon", "gauge_lat"), crs = 4326)

cirstat_points_sf <- cirstat_points_sf %>%
  mutate(
    Amax_diff_cat = cut(
      diff_days,
      breaks = c(-Inf, 14, 31, Inf),
      labels = c("0–14 days", "15–31 days", "> 31 days"),
      right  = TRUE
    )
  )

my_cols <- c(
  "0–14 days"  = "#253494",
  "15–31 days" = "#41b6c4",
  "> 31 days"  = "#edf8b1"
)

ggplot() +
  geom_sf(data = usa, fill = "grey95", color = "grey50", linewidth = 0.2) +
  geom_sf(data = states, fill = NA, color = "grey50", linewidth = 0.15) +
  geom_sf(
    data  = cirstat_points_sf,
    aes(fill = Amax_diff_cat),
    shape = 21,
    color = "grey35",
    size  = 2.8, 
    stroke = 0.15
  ) +
  scale_fill_manual(
    values = my_cols,
    name   = "Difference in mean flood occurrence date\n(Observed vs. ISIMIP)",
    drop   = FALSE
  ) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.justification = "left",
    legend.box.just = "left",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

ggsave("Figures/seas1.png", width = 8, height = 5, dpi = 600)

#Spearman correlation
Amax <- read.csv("Output/ISIMIP/ISIMIP_output/Amax_df.csv")

Amax <- Amax %>%
  mutate(id = sprintf("%08d", id))

Amax_points_df <- Amax %>%
  left_join(topo2, by = "id")

Amax_points_sf <- st_as_sf(Amax_points_df, coords = c("gauge_lon", "gauge_lat"), crs = 4326)

ggplot() +
  geom_sf(data = usa, fill = "grey95", color = "grey50", linewidth = 0.2) +
  geom_sf(data = states, fill = NA, color = "grey50", linewidth = 0.15) +
  geom_sf(data = Amax_points_sf, aes(fill = r_amax_obs_isi_sp), shape = 21, color = "grey35", size = 2.8, stroke = 0.15) +
  scale_fill_stepsn(
    colours = brewer.pal(9, "YlGnBu"),
    limits  = c(0, 1),
    breaks  = seq(0, 1, 0.2),
    name    = "Spearman correlation of flood seasonality (Observed vs. ISIMIP)",
    guide   = guide_coloursteps(
      title.position = "top",
      label.position = "bottom",
      even.steps = TRUE,
      barwidth = unit(8, "cm"),
      barheight = unit(0.5, "cm"))
  ) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.justification = "left",
    legend.box.just = "left",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

ggsave("Figures/seas2.png", width = 8, height = 5, dpi = 600)

#R (Floodconcentration)
my_blues <- c("#f7fbff", "#c6dbef","#6baed6", "#2171b5", "#08306b")

ggplot() +
  geom_sf(data = usa, fill = "grey95", color = "grey50", linewidth = 0.2) +
  geom_sf(data = states, fill = NA, color = "grey50", linewidth = 0.15) +
  geom_sf(data = cirstat_points_sf, aes(fill = R_obs), shape = 21, color = "grey40", size = 2.8, stroke = 0.15) +
  scale_fill_stepsn(
    colours = my_blues,
    limits  = c(0, 1),
    breaks  = seq(0, 1, 0.2),
    name    = "Flood timing concentration (R) (Observed)", 
    guide   = guide_coloursteps(
      title.position = "top",
      label.position = "bottom",
      even.steps = TRUE,
      barwidth = unit(8, "cm"),
      barheight = unit(0.5, "cm"))
  ) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.justification = "left",
    legend.box.just = "left",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

ggsave("Figures/seas3.png", width = 8, height = 5, dpi = 600)

##Magnitude
mag_sum <- read.csv("Output/ISIMIP/ISIMIP_output/mag_sum.csv")

mag_sum <- mag_sum %>%
  mutate(id = sprintf("%08d", id))

mag_points_df <- mag_sum %>%
  left_join(topo2, by = "id")

mag_points_sf <- st_as_sf(mag_points_df, coords = c("gauge_lon", "gauge_lat"), crs = 4326)

mag_points_sf <- mag_points_sf %>%
  mutate(
    mag_diff_cat = cut(
      diff_mag,
      breaks = c(-100, -75, -50, -25, 0, 25, 50, Inf),
      labels = c(
        "-100 to -75",
        "-75 to -50",
        "-50 to -25",
        "-25 to 0",
        "0 to 25",
        "25 to 50",
        "> 50"
      ),
      include.lowest = TRUE,
      right = TRUE
    )
  )

my_cols <- brewer.pal(9, "BrBG")

ggplot() +
  geom_sf(data = usa,    fill = "grey95", color = "grey40", linewidth = 0.3) +
  geom_sf(data = states, fill = NA,       color = "grey40", linewidth = 0.25) +
  geom_sf(
    data  = mag_points_sf,
    aes(fill = mag_diff_cat),
    shape = 21,   
    color = "grey35",
    size  = 2.8, 
    stroke = 0.15
  ) +
  scale_fill_manual(
    values = my_cols,
    name   = paste(
      "Relative difference in annual maximum flood magnitude\n(ISIMIP − Observed) [%]"),
    drop = FALSE,
    na.value = "grey60"
  )+
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.justification = "left",
    legend.box.just = "left",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

ggsave("Figures/seas4_mag.png", width = 8, height = 5, dpi = 600)

#04.2 ISIMIP climate variables bias 
clim_dev_df <- read.csv("Output/ISIMIP/ISIMIP_output/clim_dev_df.csv")

#histogram
#define breaks
bias_breaks <- seq(-60, 60, by = 5)

ggplot(clim_dev_df, aes(x = Pbias_cum_pct)) +
  geom_histogram(
    breaks = bias_breaks,
    fill   = "grey40",
    color  = "white"
  ) +
  labs(
    #title = "Distribution of cumulative precipitation bias (ISIMIP vs. CAMELS)",
    x     = "Precipitation bias [%]",
    y     = "Number of catchments"
  ) +
  scale_x_continuous(
    breaks = seq(-60, 60, by = 10),
    limits = c(-60, 60)
  )+
  coord_cartesian(xlim = c(-60, 60)) +
  theme_minimal(base_size = 12)

ggsave("Figures/hist_prec_isimip.png", width = 8, height = 5, dpi = 600)

bias_sum <- clim_dev_df %>%
  mutate(id = sprintf("%08d", as.integer(id)))

bias_points_df <- bias_sum %>%
  left_join(topo2, by = "id")

bias_points_sf <- st_as_sf(bias_points_df, coords = c("gauge_lon", "gauge_lat"), crs = 4326)

#define bias classes 
bias_points_sf <- bias_points_sf %>%
  mutate(
    Pbias_cat = cut(
      Pbias_cum_pct,
      breaks = c(-Inf, -50, -30, -20, -10, 10, 20, 30, 50),
      labels = c(
        "< -50",
        "-50 to -30",
        "-30 to -20",
        "-20 to -10",
        "-10 to 10",
        "10 to 20",
        "20 to 30",
        "30 to 50"
      ),
      include.lowest = TRUE,
      right = TRUE
    )
  )


my_cols  <- c(
  "< -50"       = "#543005",
  "-50 to -30"  = "#8c510a",
  "-30 to -20"  = "#bf812d",
  "-20 to -10"  = "#dfc27d",
  "-10 to 10"    = "#f5f5f5",
  "10 to 20"    = "#80cdc1",
  "20 to 30"    = "#35978f",
  "30 to 50"    = "#01665e"
)

names(my_cols) <- levels(bias_points_sf$Pbias_cat)

ggplot() +
  geom_sf(data = usa, fill = "grey95", color = "grey50", linewidth = 0.2) +
  geom_sf(data = states, fill = NA, color = "grey50", linewidth = 0.15) +
  geom_sf(
    data  = bias_points_sf,
    aes(fill = Pbias_cat),
    shape = 21,
    color = "grey35",
    size  = 2.8, 
    stroke = 0.15
  ) +
  scale_fill_manual(
    values   = my_cols,
    name     = "Cumulative precipitation bias\n(ISIMIP vs. CAMELS) [%]",
    drop     = FALSE,
    na.value = "grey60"
  ) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.justification = "left",
    legend.box.just = "left",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )
ggsave("Figures/ISIMIP_prec_bias_cum.png", width = 8, height = 5, dpi = 600)

#mean_annual_bias
bias_points_sf <- bias_points_sf %>%
  mutate(
    mean_annual_bias_cat = cut(
      mean_annual_bias_pct,
      breaks = c(-Inf, -50, -30, -20, -10, 10, 20, 30, 50, Inf),
      labels = c(
        "< -50",
        "-50 to -30",
        "-30 to -20",
        "-20 to -10",
        "-10 to 10",
        "10 to 20",
        "20 to 30",
        "30 to 50",
        "> 50"
      ),
      include.lowest = TRUE,
      right = TRUE
    )
  )

my_cols <- c(
  "< -50"       = "#543005",
  "-50 to -30"  = "#8c510a",
  "-30 to -20"  = "#bf812d",
  "-20 to -10"  = "#dfc27d",
  "-10 to 10"    = "#f5f5f5",
  "10 to 20"    = "#80cdc1",
  "20 to 30"    = "#35978f",
  "30 to 50"    = "#01665e", 
  "> 50" = "#003c30"
)

names(my_cols) <- levels(bias_points_sf$mean_annual_bias_cat)

ggplot() +
  geom_sf(data = usa, fill = "grey95", color = "grey50", linewidth = 0.2) +
  geom_sf(data = states, fill = NA, color = "grey50", linewidth = 0.15) +
  geom_sf(
    data  = bias_points_sf,
    aes(fill = mean_annual_bias_cat),
    shape = 21,
    color = "grey35",
    size  = 2.8, 
    stroke = 0.15
  ) +
  scale_fill_manual(
    values   = my_cols,
    name     = "Mean annual precipitation bias\n(ISIMIP vs. CAMELS) [%]",
    drop     = FALSE,
    na.value = "grey60"
  ) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.justification = "left",
    legend.box.just = "left",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )


ggsave("Figures/ISIMIP_prec_bias_mean.png", width = 8, height = 5, dpi = 600)

#mean absolute bias 
bias_points_sf <- bias_points_sf %>%
  mutate(
    mae_Pbias_annual_cat = cut(
      mae_Pbias_annual_pct,
      breaks = c(10, 20, 30, 50, 70, 100, Inf),
      labels = c(
        "10 to 20",
        "20 to 30",
        "30 to 50",
        "50 to 70",
        "70 to 100",
        "> 100"
      ),
      include.lowest = TRUE,
      right = TRUE
    )
  )

my_cols <- c(
  "10 to 20"    = "#f5f5f5",
  "20 to 30"    = "#c7eae5",
  "30 to 50"    = "#80cdc1", 
  "50 to 70" = "#35978f",
  "70 to 100" = "#01665e",
  ">100" = "#003c30"
)

names(my_cols) <- levels(bias_points_sf$mae_Pbias_annual_cat)

ggplot() +
  geom_sf(data = usa, fill = "grey95", color = "grey50", linewidth = 0.2) +
  geom_sf(data = states, fill = NA, color = "grey50", linewidth = 0.15) +
  geom_sf(
    data  = bias_points_sf,
    aes(fill = mae_Pbias_annual_cat),
    shape = 21,
    color = "grey35",
    size  = 2.8, 
    stroke = 0.15
  ) +
  scale_fill_manual(
    values   = my_cols,
    name     = "Mean absolute annual precipitation bias\n(ISIMIP vs. CAMELS) [%]",
    drop     = FALSE,
    na.value = "grey60"
  ) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.justification = "left",
    legend.box.just = "left",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

ggsave("Figures/ISIMIP_prec_bias_abs.png", width = 8, height = 5, dpi = 600)

#05 Floodprocesses (CAMELS and ISIMIP)
#load csv of both classification results 
class_results_cam1 <- read.csv("Output/CAMELS/CAMELS_output/classification_results_camels_M1.csv")
class_results_cam2 <- read.csv("Output/CAMELS/CAMELS_output/classification_results_camels_M2.csv")

class_results_isi_hist <- read.csv("Output/ISIMIP/ISIMIP_output/class_results_hist_CP_camels.csv")
class_results_isi_fut <- read.csv("Output/ISIMIP/ISIMIP_output/class_results_future_2041_2070_ssp370_CP_camels.csv")
class_results_isi_fut_far <- read.csv("Output/ISIMIP/ISIMIP_output/class_results_future_2071_2100_ssp370_CP_camels.csv")

#process contribution map (spatial distribution)
class_results_sum <- class_results_cam1 %>% #change to right input class results csv
  mutate(id = sprintf("%08d", id))

points_process <- class_results_sum %>%
  left_join(topo2, by = "id")

points_process <- st_as_sf(points_process,
                               coords = c("gauge_lon", "gauge_lat"),
                               crs = 4326)
# Define labels and plotting order
process_map <- c(
  perc_rainfall     = "Short rainfall",
  perc_longrainfall = "Long rainfall",
  perc_rainandsnow  = "Rainfall/Snowmelt",
  perc_snowmelt     = "Snowmelt",
  perc_soilsat      = "Excess rainfall"
)
process_levels <- c(
  "Excess rainfall",
  "Short rainfall",
  "Long rainfall",
  "Rainfall/Snowmelt",
  "Snowmelt"
)

#Bounding box processes 
xlim <- c(-126, -66)
ylim <- c(24, 50)

#Convert process contribution columns to long format and assign process names
pts_long <- points_process %>%
  select(all_of(c("gauge_id", names(process_map)))) %>%
  pivot_longer(
    cols      = all_of(names(process_map)),
    names_to  = "process",
    values_to = "share"
  ) %>%
  mutate(
    process = factor(process_map[process], levels = process_levels)
  )

my_blues    <- colorRampPalette(brewer.pal(9, "Blues"))(20)
breaks_full <- seq(0, 100, 5)              # small boxes every 5%
labels_full <- ifelse(breaks_full %% 25 == 0, breaks_full, "")   # label ticks at 0,25,50,75,100

ggplot() +
  geom_sf(data = usa, fill = "grey95", color = "grey40", linewidth = 0.3) +
  geom_sf(data = states, fill = NA, color = "grey40", linewidth = 0.25) +
  geom_sf(
    data  = pts_long,
    aes(color = share),
    shape = 16,
    size  = 3
  ) +
  scale_color_stepsn(
    colours = my_blues,
    limits  = c(0, 100),
    breaks  = breaks_full,         
    labels  = labels_full,      
    name    = "CAMELS: Process contribution [%]", #adjust the right name 
    guide   = guide_colorbar(
      barwidth      = unit(4.5, "cm"),
      barheight     = unit(0.5, "cm"),
      title.position = "top",
      label.position = "bottom"
    )
  ) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  facet_wrap(~ process, ncol = 3) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position  = "bottom",
    legend.box       = "vertical",
    legend.title     = element_text(hjust = 0.5),
    legend.key.width = unit(2, "cm"),
    strip.text       = element_text(face = "bold")
  )

ggsave("Figures/CAMELS_flood_procM1.png", width = 12, height = 6, dpi = 300)

#scatterplot processes CAMELS vs ISIMIP
camels_df <- class_results_cam2 %>%
  mutate(id = sprintf("%08d", as.integer(id)))  

isimip_df <- class_results_isi_hist %>%
  mutate(id = sprintf("%08d", as.integer(id)))

df_join <- camels_df %>%
  inner_join(isimip_df, by = "id", suffix = c("_camels", "_isimip"))

proc_map <- c(
  perc_soilsat      = "Excess rainfall",
  perc_rainfall     = "Short rainfall",
  perc_longrainfall = "Long rainfall",
  perc_rainandsnow  = "Rainfall/Snowmelt",
  perc_snowmelt     = "Snowmelt",
  perc_noclass = "NA"
)

proc_levels <- c("Excess rainfall", "Short rainfall", "Long rainfall", "Rainfall/Snowmelt","Snowmelt")

df_long <- df_join %>%
  pivot_longer(
    cols = matches("^perc_.*_(camels|isimip)$"),
    names_to = c("process", "dataset"),
    names_pattern = "(perc_.*)_(camels|isimip)$",
    values_to = "share"
  )

df_plot <- df_join %>%
  pivot_longer(
    cols = matches("^perc_.*_(camels|isimip)$"),
    names_to = c("process", "dataset"),
    names_pattern = "(perc_.*)_(camels|isimip)$",
    values_to = "share"
  ) %>%
  mutate(
    process = recode(process, !!!proc_map),
    process = factor(process, levels = proc_levels)
  ) %>%
  filter(process != "NA") %>%
  pivot_wider(
    names_from = dataset,
    values_from = share
  )

ggplot(df_plot, aes(x = camels, y = isimip)) +
  geom_point(shape = 1, colour = "blue", alpha = 0.6, size = 1.4, stroke = 0.7) + 
  geom_abline(slope = 1, intercept = 0, linewidth = 0.4, colour = "grey40") +
  coord_equal(xlim = c(0, 100), ylim = c(0, 100)) +
  facet_wrap(~ process, nrow = 1) + 
  scale_x_continuous(
    breaks = c(0, 25, 50, 75, 100),
    minor_breaks = NULL
  ) +
  scale_y_continuous(
    breaks = c(0, 25, 50, 75, 100),
    minor_breaks = NULL
  ) +
  labs(
    x = "CAMELS process contribution [%]",
    y = "ISIMIP process contribution [%]"
  ) +
  theme_minimal(base_size = 9)+
  theme(
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    strip.text   = element_text(size = 9, face = "plain"),
    panel.grid.minor = element_blank(),               
    panel.grid.major = element_line(linewidth = 0.5),
  )

ggsave("Figures/cam_vs_isi_scatter.png", width = 12, height = 3, dpi = 600)

#Plot process contribution changes over time (ISIMIP)
hist_long <- class_results_isi_hist %>%
  select(id, starts_with("perc_")) %>%
  pivot_longer(
    cols = starts_with("perc_"),
    names_to = "process",
    values_to = "share"
  ) %>%
  mutate(
    process = sub("^perc_", "", process),
    time = "Historic")

fut_mid_long <- class_results_isi_fut %>%
  select(id, starts_with("perc_")) %>%
  pivot_longer(
    cols = starts_with("perc_"),
    names_to = "process",
    values_to = "share"
  ) %>%
  mutate(
    process = sub("^perc_", "", process),
    time = "Mid future"
  )

fut_far_long <- class_results_isi_fut_far %>%
  select(id, starts_with("perc_")) %>%
  pivot_longer(
    cols = starts_with("perc_"),
    names_to = "process",
    values_to = "share"
  ) %>%
  mutate(
    process = sub("^perc_", "", process),
    time = "Far future"
  )

plot_time <- bind_rows(
  hist_long,
  fut_mid_long,
  fut_far_long
) %>%  
  mutate(
    time = factor(time,
                  levels = c("Historic","Mid future","Far future"), labels = c("Historic", "Near future", "Future"))
  )

process_order <- c(
  "soilsat",       
  "rainfall",       
  "longrainfall",   
  "rainandsnow",  
  "snowmelt",
  "noclass"
)

plot_time <- plot_time %>%
  mutate(
    process = factor(
      process,
      levels = process_order,
      labels = c(
        "Excess rainfall",
        "Short rainfall",
        "Long rainfall",
        "Rainfall / Snowmelt",
        "Snowmelt",
        "No class"
      )
    )
  )

plot_time_active <- plot_time %>%
  group_by(id, process) %>%
  filter(any(share > 0)) %>%
  ungroup()

ggplot(plot_time_active, aes(x = time, y = share, group = id)) +
  geom_line(alpha = 0.08, linewidth = 0.3) +
  geom_point(alpha = 0.15, size = 0.4) +
  stat_summary(
    aes(group = 1),
    fun = median, geom = "line",
    linewidth = 0.4, color = "#D55E00"
  ) +
  stat_summary(
    aes(group = 1),
    fun = median, geom = "point",
    size = 0.6, color = "#D55E00"
  ) +
  facet_wrap(~ process) +
  labs(x = "", y = "Process contribution (%)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 8),
    strip.text  = element_text(size = 9)
  )

ggsave("Figures/proc_con_ssp370.png", width = 12, height = 6, dpi = 600)

#Table: dominant process change historical compared to future ISIMIP
flood_proc_hist <- read.csv("Output/ISIMIP/ISIMIP_output/flood_proc_hist_CP_camels.csv")
flood_proc_fut <- read.csv("Output/ISIMIP/ISIMIP_output/flood_proc_future_2071_2100_ssp585_CP_camels.csv")

process_change <- flood_proc_hist %>%
  left_join(flood_proc_fut, by = "id")

process_trans <- process_change %>%
  count(Dom_flood_proc_isi.x, Dom_flood_proc_isi.y)

process_trans
