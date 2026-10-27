# Example validation plot: Observed vs AntTerrClim
set.seed(42)

# Simulated example data (replace with your own)
observed <- runif(1000, min = -25, max = 5)           # e.g. temperature in °C
antterrclim <- observed + rnorm(1000, mean = 0, sd = 2)  # predictions with some error

# Plot
plot(observed, antterrclim,
     pch  = 16,                 # solid dots
     cex  = 0.5,                # point size
     col  = "black",
     xlab = "Observed (°C)",
     ylab = "AntTerrClim (°C)",
     main = "Temperature",
     bty  = "l")                # L-shaped box (axes on left and bottom only)

# Dashed 1:1 line
abline(a = 0, b = 1, lty = 2, col = "black")



# Example validation plot: Observed vs AntTerrClim, with some bias
set.seed(42)

# Simulated example data (replace with your own)
observed <- runif(1000, min = -25, max = 5)   # e.g. temperature in °C

bias_offset <- 1.5    # constant warm bias (°C)
bias_slope  <- 0.9    # <1 means the model under-represents the extremes
noise_sd    <- 2      # random scatter

antterrclim <- bias_offset + bias_slope * observed + rnorm(1000, 0, noise_sd)

# Save as PNG
png("antterrclim_validation.png",
    width = 1600, height = 1600, res = 300)   # 1600 px at 300 dpi, i.e. ~5.3 x 5.3 inches

plot(observed, antterrclim,
     pch  = 16,
     cex  = 0.5,
     col  = "black",
     xlab = "Observed (°C)",
     ylab = "AntTerrClim (°C)",
     main = "Temperature",
     bty  = "l")

# Dashed 1:1 line
abline(a = 0, b = 1, lty = 2, col = "black")

dev.off()   # closes the file and writes the PNG



# Example validation plot: Observed vs AntTerrClim, sized for a small figure
set.seed(42)

# Simulated example data (replace with your own)
n <- 600                                       # fewer points so they don't clump when small
observed <- runif(n, min = -25, max = 5)       # e.g. temperature in °C

bias_offset <- 1.5
bias_slope  <- 0.9
noise_sd    <- 2

antterrclim <- bias_offset + bias_slope * observed + rnorm(n, 0, noise_sd)

# Save as PNG at the actual size it will appear on the page
png("antterrclim_validation.png",
    width = 2.5, height = 2.5, units = "in",   # set to the final printed size
    res = 600,                                 # high res so it stays sharp
    pointsize = 8)                             # base text size in points

# Tight margins (bottom, left, top, right) and closer axis labels
par(mar = c(3, 3, 1.5, 0.5),
    mgp = c(1.8, 0.6, 0))

plot(observed, antterrclim,
     pch      = 16,
     cex      = 0.4,          # small points to limit overplotting
     col      = "grey50",
     xlab     = "Observed (°C)",
     ylab     = "AntTerrClim (°C)",
     cex.lab  = 1,            # axis label size
     cex.axis = 0.9,          # tick label size
     font.main = 1,           # plain (not bold) title
     bty      = "l",
     lwd      = 1)            # axis line weight

# Dashed 1:1 line, a bit thicker so it survives shrinking
abline(a = 0, b = 1, lty = 2, lwd = 1.8, col = "black")

dev.off()

