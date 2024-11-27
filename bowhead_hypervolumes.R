library(data.table)
library(doSNOW)
library(pbapply)
library(hypervolume)
library(terra)
library(raster)
library(sf)
library(gtools)

devtools::source_gist("https://gist.github.com/scbrown86/87a8c49117d0be7ea8730a6f4d17fa2c")

#### FULL HYPERVOLUME ####
# Read in the full niche (fossils matched to climate/env records)
bm_cal <- setDT(readRDS("Data/Fossils/bowhead_hypervolume_fossil_records_matched.RDS")[[1]])
bm_cal
keep_cols <- colnames(bm_cal)[c(1:8, 10,13,16)]
bm_cal <- bm_cal[, ..keep_cols]
bm_cal

bm_val <- setDT(readRDS("Data/Fossils/bowhead_hypervolume_fossil_records_matched.RDS")[[2]])
bm_val <- bm_val[, ..keep_cols]
bm_val

# Define the climate and ID columns
## exclude Slope and Rugosity due to lots of NA
climvars <- colnames(bm_cal)[c(9:11)]
id_cols <- colnames(bm_cal)[c(1:8)]

# PCA
bm_full <- rbindlist(list(bm_cal, bm_val))
env_pca <- prcomp(as.matrix(bm_full[, ..climvars]), center = TRUE,
                  scale. = TRUE, retx = TRUE)
env_pca; summary(env_pca)
cols.group <- character(nrow(bm_full))
cols.group[] <- "#1b9e77"
cols.group[bm_full$Region == "Svalbard"] <- "#d95f02"
pch.group <- numeric(nrow(bm_full))
pch.group[] <- 1
# plot the PCA
# png(filename = "Figs/pca_bowhead_summer.png",
#     width = 8, height = 8, units = "in",
#     type = "cairo-png",
#     res = 320, antialias = "subpixel")
{plot(env_pca$x[,1], env_pca$x[,2],
     xlab=paste("PCA 1 (", round(summary(env_pca)$importance[2]*100, 1), "%)", sep = ""),
     ylab=paste("PCA 2 (", round(summary(env_pca)$importance[5]*100, 1), "%)", sep = ""),
     xlim = c(-6,6), ylim = c(-6, 6),
     # pch = pch.group,
     col = cols.group,
     bg = cols.group, cex = 1.5, las = 1)
# overlay validation points
points(env_pca$x[bm_full$Use == "Validation", 1],
       env_pca$x[bm_full$Use == "Validation", 2],
       col = "black", bg = "#7570b3", pch = 23, cex = 1.5)
abline(v = 0, lty = 2, col = "grey50")
abline(h = 0, lty = 2, col = "grey50")
# Get co-ordinates of variables (loadings), and multiply by 5
l.x <- env_pca$rotation[, 1]*5
l.y <- env_pca$rotation[, 2]*5
arrows(x0 = 0, x1 = l.x, y0 = 0, y1 = l.y, col = "black", length = 0.15, lwd = 3)
# Label position
l.pos <- l.y # Create a vector of y axis coordinates
lo <- which(l.y < 0) # Get the variables on the bottom half of the plot
hi <- which(l.y > 0) # Get variables on the top half
# Replace values in the vector
l.pos <- replace(l.pos, lo, "1")
l.pos <- replace(l.pos, hi, "3")
# Variable labels
l.pos[1] <- 2
l.pos[2] <- 2
l.pos[3] <- 4
text(l.x, l.y,
     #labels = row.names(env_pca$rotation),
      labels = c("SST", "Sal.", "SIC"),
     #col = "#1f78b4",
     col = "black",
     pos = l.pos, cex = 1.5)
# Get individuals (observations) as a matrix
tab <- matrix(c(env_pca$x[,1], env_pca$x[,2]), ncol=2)
# Calculate correlations
c1 <- cor(tab[bm_full$Region == "Svalbard",], method = "kendall")
c2 <- cor(tab[bm_full$Region != "Svalbard",], method = "kendall")
polygon(ellipse::ellipse(c1*(max(abs(env_pca$rotation))*1),
                centre = colMeans(tab[bm_full$Region == "Svalbard",]), level = 0.95),
        col = adjustcolor("#d95f02", alpha.f = 0.25), border = "#d95f02")
polygon(ellipse::ellipse(c2*(max(abs(env_pca$rotation))*1),
                centre = colMeans(tab[bm_full$Region != "Svalbard",]), level = 0.95),
        col = adjustcolor("#1b9e77", alpha.f = 0.25), border = "#1b9e77")
legend("bottomright",
       legend = c("Canada", "Svalbard", "validation"),
       col = c("#1b9e77", "#d95f02", "black"), pt.bg = c("#1b9e77", "#d95f02", "#7570b3"),
       pch = c(1, 1, 23), pt.cex = 1.5)
}
# dev.off()

# quick plots
# png(filename = "Figs/pairs_plot.png",
#     width = 8, height = 8, units = "in",
#     type = "cairo-png",
#     res = 320, antialias = "subpixel")
pairs(data.frame(SST = bm_full[["SST_summ"]],
                 SIC = bm_full[["IceCon_summ"]],
                 Sal. = bm_full[["Sal_summ"]]),
      col = cols.group, pch = ifelse(bm_full[["Use"]] == "Calibration", 0, 24))
# dev.off()
# center and scale measurements
hv_vars <- climvars
df_z <- scale(bm_full[, ..hv_vars], center = TRUE, scale = TRUE)
tmp <- cor(df_z, method = "kendall", use = "complete.obs")
tmp[!lower.tri(tmp)] <- NA
print(tmp, digits = 3, na.print = "")
#               SST_summ Sal_summ IceCon_summ
# SST_summ
# Sal_summ       0.272
# IceCon_summ   -0.489   -0.193

summary(bm_full[, ..hv_vars])

# Quick test
hv_list <- new("HypervolumeList")
hv_list@HVList = vector(mode = "list", length = 4)
hv_list@HVList[[1]] <- hypervolume_gaussian(data = df_z[bm_full[["Use"]] == "Calibration", ],
                                            name = "Calibration", samples.per.point = 3,
                                            kde.bandwidth = estimate_bandwidth(df_z[bm_full[["Use"]] == "Calibration", ], method = "silverman-1d"),
                                            verbose = TRUE)
hv_list@HVList[[2]] <- hypervolume_gaussian(data = df_z[bm_full[["Use"]] == "Validation", ],
                                            name = "Validation", samples.per.point = 3,
                                            kde.bandwidth = estimate_bandwidth(df_z[bm_full[["Use"]] == "Validation", ], method = "silverman-1d"),
                                            verbose = TRUE)
hv_list@HVList[[3]] <- hypervolume_gaussian(data = df_z[bm_full[["Region"]] == "Canada", ],
                                            name = "Canada", samples.per.point = 3,
                                            kde.bandwidth = estimate_bandwidth(df_z[bm_full[["Region"]] == "Canada", ], method = "silverman-1d"),
                                            verbose = TRUE)
hv_list@HVList[[4]] <- hypervolume_gaussian(data = df_z[bm_full[["Region"]] == "Svalbard", ],
                                            name = "Svalbard", samples.per.point = 3,
                                            kde.bandwidth = estimate_bandwidth(df_z[bm_full[["Region"]] == "Svalbard", ], method = "silverman-1d"),
                                            verbose = TRUE)
saveRDS(hv_list, "Outputs/test_hypervolumes.RDS", compress = TRUE)

plot(hv_list,
     num.points.max.data = 3000, num.points.max.random = 4000)

plot(hv_list[[c(3:4)]], cex.data = 1,
     col = c("Orange", "Green"),
     num.points.max.data = 3000, num.points.max.random = 4000)

# Make a test projection of environmental suitability
clim_dir <- "path/to/climate_data/"
ras_ssts <- crop(brick(file.path(clim_dir, "bowhead_env_vars_11700BP_2020CE.nc"), varname = "tos")[[11737]],
                 c(-180,180,50,90))
ras_icec_s <- crop(brick(file.path(clim_dir, "bowhead_env_vars_11700BP_2020CE.nc"), varname = "si")[[11737]],
                   c(-180,180,50,90))
ras_sals <- crop(brick(file.path(clim_dir, "bowhead_env_vars_11700BP_2020CE.nc"), varname = "so")[[11737]],
                 c(-180,180,50,90))
s <- stack(ras_ssts, ras_sals, ras_icec_s)
names(s) <- colnames(hv_list@HVList[[1]]@Data)
s <- scale(s, center = attr(df_z,"scaled:center"), scale = attr(df_z,"scaled:scale")); s
plot(s, nc = 1)

# projections using default settings
cal_prj <- hypervolume_project(hv_list@HVList[[1]], rasters = s, parallel = TRUE, n.cores = 3)
val_prj <- hypervolume_project(hv_list@HVList[[2]], rasters = s, parallel = TRUE, n.cores = 3)
can_prj <- hypervolume_project(hv_list@HVList[[3]], rasters = s, parallel = TRUE, n.cores = 3)
sval_prj <- hypervolume_project(hv_list@HVList[[4]], rasters = s, parallel = TRUE, n.cores = 3)
prj_stack <- stack(cal_prj, val_prj, can_prj, sval_prj)
names(prj_stack) <- c("Calibration", "Validation", "Canada", "Sval")
spplot(prj_stack)

# bandwidth
bw_tune <- estimate_bandwidth(df_z[bm_full[["Use"]] == "Calibration", ],
                              method = "cross-validation")
bw_tune

# grid for tuning the hypervolume settings
tune_grid <- expand.grid(
  quant = c(0.95, 0.975, 0.99),
  kde_mul = seq(1.0, 1.5, 0.25),
  sd.count = c(3, 4))
tune_grid

# HypervolumeList class for storing the results
hv_list <- new("HypervolumeList")
hv_list@HVList = vector(mode = "list")

# Process the tuning in parallel
cls <- makeCluster(9L)
registerDoSNOW(cls)
pb <- txtProgressBar(max = nrow(tune_grid), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Iterate through the hypervolume settings
hv_list@HVList <- foreach(row = seq_len(nrow(tune_grid)),
                          .packages = c("hypervolume", "data.table"),
                          .inorder = TRUE, .verbose = FALSE,
                          .multicombine = FALSE,
                          .errorhandling = "remove",
                          .noexport = ls(globalenv()),
                          .export = c("tune_grid", "df_z", "bw_tune", "bm_full"),
                          .options.snow = opts) %dopar% {
                            qt <- tune_grid[row, 1]
                            mul <- tune_grid[row, 2]
                            sdc <- tune_grid[row, 3]
                            hv <- hypervolume(df_z[bm_full[["Use"]] == "Calibration", ],
                                              method = "gaussian",
                                              name = sprintf("summer_hv_q_%s_bw_%s_sd_%s",
                                                             as.character(qt),
                                                             as.character(mul),
                                                             as.character(sdc)),
                                              chunk.size = round(nrow(df_z[bm_full[["Use"]] == "Calibration", ])/6),
                                              verbose = FALSE,
                                              samples.per.point = 5,
                                              quantile.requested = qt,
                                              kde.bandwidth = estimate_bandwidth(df_z, method = "fixed", bw_tune*mul),
                                              sd.count = sdc)
                            return(hv)
                          }
registerDoSEQ(); stopCluster(cls)
summary(hv_list)
saveRDS(hv_list, "Outputs/summer_3d_model_hv_bandwidth_tuning.RDS", compress = TRUE)

# png(filename = "Figs/hypervolume_bandwidth_estimates_summer3d.png",
#     width = 8, height = 8, units = "in", res = 320,
#     type = "cairo")
plot(hv_list,
     show.centroid = TRUE, show.random = TRUE, show.data = TRUE,
     contour.type = "kde",
     point.alpha.min = 0.0001,
     cex.legend = 0.8,
     contour.lwd = 1,
     plot.function.additional = function(i, j) {
       points(x = df_z[bm_full[["Use"]] != "Calibration", i],
              y = df_z[bm_full[["Use"]] != "Calibration", j], col = "black")
     })
# dev.off()

# subset hypervolumes to 90th percentile of volumes
## selects the biggest hypervolumes
q90 <- quantile(get_volume(hv_list), 0.90)
idx <- which(lapply(hv_list@HVList, function(x) get_volume(x) >= q90) == TRUE)
subset_hv <- hv_list[[idx]]
names(subset_hv@HVList) <- sapply(subset_hv@HVList, function(x) x@Name)
subset_hv

get_volume(subset_hv)
# summer_hv_q_0.99_bw_1.5_sd_3 summer_hv_q_0.975_bw_1.5_sd_4
#                     9.212611                      9.294653

# % of validation points included in the hypervolume
points_in <- sapply(lapply(subset_hv@HVList, function(x) {
  tot <- hypervolume_inclusion_test(x, df_z[bm_full[["Use"]] == "Validation", ],
                                    fast.or.accurate = "accurate",
                                    verbose = FALSE)
  per <- round(sum(tot)/length(tot),2)
  setattr(x, name = "val_per", value = per)
  return(per)
}), "[", 1)
points_in

# summer_hv_q_0.99_bw_1.5_sd_3 summer_hv_q_0.975_bw_1.5_sd_4
#                         0.97                          0.96

# Plot of hypervolumes with validation points
# png(filename = "Figs/hypervolume_bandwidth_estimates_summer3d_wValPoints.png",
#     width = 8, height = 8, units = "in", res = 320,
#     type = "cairo")
plot(subset_hv,
     colors = c("skyblue2", "gold"),
     show.centroid = TRUE, show.random = TRUE, show.data = TRUE,
     num.points.max.random = 5000,
     num.points.max.data = nrow(df_z),
     point.dark.factor = 0.5,
     contour.type = "kde", reshuffle = TRUE,
     plot.function.additional = function(i,j) {
       points(x = df_z[bm_full[["Use"]] == "Validation", i],
              y = df_z[bm_full[["Use"]] == "Validation", j],
              col = "black", cex = 1, pch = 16)
     })
# dev.off()

# combinations for overlap
## we only selected 2 hypervolumes, but this code will work if we selected n hypervolumes
p_names <- 1:length(subset_hv@HVList)
combos <- combinations(n = length(p_names), r = 2, v = p_names, set = TRUE, repeats.allowed = FALSE)

overlap <- matrix(NA, nrow = nrow(combos), ncol = 3)
overlap[, 1] <- combos[, 1]
overlap[, 2] <- combos[, 2]

for (i in 1:nrow(overlap)) {
  s1 <- overlap[i, 1]; s2 <- overlap[i, 2]
  this_set <- hypervolume_set(hv_list@HVList[[s1]], hv_list@HVList[[s2]], check.memory = FALSE,
                             verbose = FALSE)
  # calculate a Sorensen overlap index (2 x shared volume / sum of |hv1|  |hv2|)
  overlap[i, 3] <- hypervolume_overlap_statistics(this_set)["sorensen"]
}
overlap # ~82% overlap

hypervolume_distance(subset_hv@HVList[[1]], subset_hv@HVList[[2]],
                     type = "centroid",
                     num.points.max = 20000, check.memory = FALSE)

attr(subset_hv@HVList[[1]], "val_per"); attr(subset_hv@HVList[[2]], "val_per")

## summer_hv_q_0.99_bw_1.5_sd_3, biggest HV with most number of validation points ###

#### PLOTS AND TUNING PROJECTIONS ####

full_hv <- subset_hv@HVList[["summer_hv_q_0.99_bw_1.5_sd_3"]]
full_hv
summary(full_hv); plot(full_hv, show.legend = FALSE, contour.type = "kde")
setattr(full_hv@Data, "scaled:center", attr(df_z, "scaled:center"))
setattr(full_hv@Data, "scaled:scale", attr(df_z, "scaled:scale"))
attr(full_hv@Data, "scaled:center")
attr(full_hv@Data, "scaled:scale")
colnames(full_hv@Data)
saveRDS(full_hv, "Outputs/full_hv_summermodel_3d.RDS", compress = TRUE)

# Make a test projection of environmental suitability
ras_ssts <- brick(file.path(clim_dir, "bowhead_env_vars_11700BP_2020CE.nc"), varname = "tos")
ras_icec_s <- brick(file.path(clim_dir, "bowhead_env_vars_11700BP_2020CE.nc"), varname = "si")
ras_sals <- brick(file.path(clim_dir, "bowhead_env_vars_11700BP_2020CE.nc"), varname = "so")

# contemporary range
bm_range <- st_read("Data/Shapefiles/bowhead_range_iucn.shp")
bm_range <- bm_range[bm_range$BINOMIAL == "Balaena mysticetus", ]
bm_range
plot(bm_range[1])

s <- readAll(stack(ras_ssts[["X.55"]], ras_sals[["X.55"]], ras_icec_s[["X.55"]]))
s <- crop(s, c(-180,180,50,90))
names(s) <- colnames(full_hv@Data)
plot(s, addfun = function() lines(as_Spatial(bm_range)))

# scale relative to full niche
s <- scale(s,
           center = attr(full_hv@Data, "scaled:center"),
           scale = attr(full_hv@Data,"scaled:scale"))
s; spplot(s)

tst_prj <- hypervolume_project(full_hv, rasters = s,
                               verbose = FALSE, parallel = TRUE,
                               n.cores = 3L, edges.zero.distance.factor = 3,
                               weight.exponent = -1, set.edges.zero = TRUE)
spplot(tst_prj)

def_prj <- hypervolume_project(full_hv, rasters = s,
                               parallel = TRUE,
                               n.cores = 3L,
                               verbose = FALSE,
                               set.edges.zero = TRUE)
def_prj
spplot(def_prj)

spplot(stack(def_prj, tst_prj))
plot(tst_prj)
points(bm_full[Use == "Validation", c("Lon", "Lat")], pch = 1, col = "red")
points(20.5,80.5, pch = 16, col = "black")
extract(tst_prj, bm_full[Use == "Validation", c("Lon", "Lat")], method = "bilinear", na.rm = TRUE)

# Tuning grid for projections
tune_grid <- expand.grid(
  edges = seq(1, 5, 0.5),
  weight = -seq(1, 3, 1))
tune_grid

projection_tuning <- pbsapply(1:nrow(tune_grid), function(i) {
  # probabilities
  cal_hs <- hypervolume_estimate_probability(hv = full_hv,
                                             points = df_z[bm_full[["Use"]] == "Calibration", ],
                                             parallel = TRUE, n.cores = 3L,
                                             weight.exponent = tune_grid[i, 2],
                                             set.edges.zero = TRUE,
                                             edges.zero.distance.factor = tune_grid[i, 1],
                                             verbose = FALSE) # calibration points
  val_hs <- hypervolume_estimate_probability(full_hv,
                                             points = df_z[bm_full[["Use"]] == "Validation", ],
                                             parallel = TRUE, n.cores = 3L,
                                             weight.exponent = tune_grid[i, 2],
                                             set.edges.zero = TRUE,
                                             edges.zero.distance.factor = tune_grid[i, 1],
                                             verbose = FALSE) # validation points
  # extract values from random points (equivalent to background points)
  ## see ?`Hypervolume-class` for details
  ## 10,000 or total (smaller or two) of background points sampled
  bg_scores <- background_probability(full_hv, n.points = 10000,
                                      random.seed = 20220510,
                                      weight.exponent = tune_grid[i, 2],
                                      set.edges.zero = TRUE,
                                      edges.zero.distance.factor = tune_grid[i, 1],
                                      verbose = FALSE)
  invisible({print({par(mfrow = c(2,2))
    hist(cal_hs, breaks = 20, xlim = range(cal_hs))
    abline(v = quantile(cal_hs, probs = c(0.10, 0.90)), col = "red", lty = 2)
    hist(val_hs, breaks = 20, xlim = range(cal_hs))
    abline(v = quantile(cal_hs, probs = c(0.10, 0.90)), col = "red", lty = 2)
    hist(bg_scores, breaks = 20, xlim = range(cal_hs))
    abline(v = quantile(cal_hs, probs = c(0.10, 0.90)), col = "red", lty = 2)
    plot(0,0, xaxt = "n", yaxt = "n",
         pch = NA, bty = "n",
         main = paste("weight = ",tune_grid[i, 2]," dist.fac = ",tune_grid[i, 1]))
    par(mfrow = c(1,1))})})
  # p10 and p95 thresholds
  p10 <- quantile(cal_hs, 0.10)
  p95 <- quantile(cal_hs, 0.95)
  # Threshold. Anything <= p10 becomes NA, anything > p95 becomes p95
  ## thresholded scores are then rescaled to {0, 1}
  cal_hs[cal_hs <= p10] <- NA; cal_hs[cal_hs > p95] <- p95; cal_hs <- scales::rescale(cal_hs, from = c(p10, p95), to = c(0, 1))
  val_hs[val_hs <= p10] <- NA; val_hs[val_hs > p95] <- p95; val_hs <- scales::rescale(val_hs, from = c(p10, p95), to = c(0, 1))
  bg_scores[bg_scores <= p10] <- NA; bg_scores[bg_scores > p95] <- p95; bg_scores <- scales::rescale(bg_scores, from = c(p10, p95), to = c(0, 1))
  invisible({print({par(mfrow = c(2,2))
    hist(cal_hs, breaks = 20, xlim = range(cal_hs, na.rm = TRUE))
    hist(val_hs, breaks = 20, xlim = range(cal_hs, na.rm = TRUE))
    hist(bg_scores, breaks = 20, xlim = range(cal_hs, na.rm = TRUE))
    plot(0,0, xaxt = "n", yaxt = "n",
         pch = NA, bty = "n",
         main = paste("weight = ",tune_grid[i, 2]," dist.fac = ",tune_grid[i, 1]))
    par(mfrow = c(1,1))})})
  ## scores
  res_list <- list(
    calAUC = enmSdmX::evalAUC(cal_hs, bg_scores, na.rm = TRUE,
                              contrastWeight = rep(1 / length(bg_scores), times = length(bg_scores))),
    valAUC = enmSdmX::evalAUC(val_hs, bg_scores, na.rm = TRUE,
                              contrastWeight = rep(1 / length(bg_scores), times = length(bg_scores))),
    # Boyce
    calBoyce = enmSdmX::evalContBoyce(cal_hs, bg_scores, na.rm = TRUE, binWidth = 0.1,
                                      contrastWeight = rep(1 / length(bg_scores), times = length(bg_scores))),
    valBoyce = enmSdmX::evalContBoyce(val_hs, bg_scores, na.rm = TRUE, binWidth = 0.1,
                                      contrastWeight = rep(1 / length(bg_scores), times = length(bg_scores)))
  )
  return(res_list)
})
projection_tuning[] <- lapply(projection_tuning,round,2)
colnames(projection_tuning) <- paste0("E_", tune_grid[, 1],"_W_", tune_grid[,2])
projection_tuning[c(2,4), ]

# Project each layer as raster
cls <- makeCluster(3L)
clusterExport(cls, c("full_hv", "ras_ssts", "ras_icec_s", "ras_sals"))

idx <- which(getZ(ras_ssts) >= -11700)

prj_list <- pblapply(idx, function(slice, ...) {
  require(raster); require(hypervolume)
  yr <- getZ(ras_ssts)[slice]
  s <- readAll(stack(ras_ssts[[slice]], ras_sals[[slice]], ras_icec_s[[slice]]))
  names(s) <- colnames(full_hv@Data)
  ## Scale and centre values relative to full niche!
  s <- scale(s,
             center = attr(full_hv@Data,"scaled:center"),
             scale = attr(full_hv@Data,"scaled:scale"))
  prj <- hypervolume_project(full_hv,
                             rasters = s,
                             verbose = FALSE,
                             parallel = TRUE, n.cores = 3L,
                             ## from tuning
                             weight.exponent = -3, edges.zero.distance.factor = 3,
                             set.edges.zero = TRUE)
  names(prj) <- paste0("X", yr)
  writeRaster(prj,
              sprintf("Outputs/summer_model_projected_default_%s.grd", names(prj)),
              overwrite = TRUE)
  return(prj)
}, cl = cls)
stopCluster(cls)
prj_list
