library(wholebrain)
library(rgl)
library(magick)
library(EBImage)

# Set this system variable to avoid RStudio from crashing
Sys.setenv(LIBGL_ALWAYS_SOFTWARE=1)
setwd("~/FFPE/analysis/doc/")
source("../doc/glassbrain_custom.R")

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## ------------------------------------------------------------------------------------------------------
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# The code beow was used to generate the registration object but you
# can skip these steps and load the finished registration (see below)

imhe <- "../data/registration/images/FF_mouse_brain_edited.tif"

# this is the filter we used it is a list object with the following parameters.
myfilter <- list(alim = c(1, 50), #area limits for what to consider as cell bodies c(minimum, maximum)
               threshold.range = c(10000, 60000), #threshold range in the minimum and maximum in fluorescent intensity where the algorithm should search for neurons.
               eccentricity = 1000, #eccentricity (elongation) of contours sets how round you want cell bodies to be. Default is 1000 and smaller values equal to more round.
               Max = 60000, #Maximum value to display in the 8-bit rendered (sets sort of brightness contrast)
               Min = 0, #Minimum value to display in the 8-bit rendered (sets sort of brightness contrast)
               brain.threshold = 2800, #the exact value where you want to start segmeting the brain outline in autofluorescence.
               resize = 0.2, #0.28, resize parameter to match the atlas to your pixel resolution, should be between 0.03 and 0.2 for most applications.
               blur = 13, #blur parameter that sets the smoothness of the tissue outline, if magaded or jagged edges increase. Using a value fo 4 is usually recommended.
               downsample = 0.5 #downsample, default is set to 0.25 and images with a size of 15000 x 8000 pixels can then usually be run smoothly
)

# segmentation
quartz()
seg <- segment(imhe, filter = myfilter, display = TRUE)

# get the fluorescent intensity of each hematoxlyin stained nuclei as a gray-scale value
fluorescent.intensity <- scales::rescale(seg$soma$intensity)

# plot the segmentation
plot(seg$soma, pch = 16, cex = 0.25, asp = 1,
     ylim = rev(range(seg$soma$y)), col = gray(fluorescent.intensity),
     axes = F, ylab = '', xlab = '')


# run first pass at rgeistration
quartz()
coord <- -2.2
regi <- registration(imhe, filter = myfilter, coordinate = coord, right.hemisphere = TRUE)

# Add correspondance points (iterative this step until the registration is satisfactory)
regi <- add.corrpoints(regi)

# update the registration
quartz()
regi <- registration(imhe, filter = myfilter, coordinate = coord, correspondance = regi, right.hemisphere = TRUE)


## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## ------------------------------------------------------------------------------------------------------
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Load the prepared registration object with manually added correspondance points
regi <- readRDS("../data/registration/registration_HE_FFPE_mouse_brain")

# Extract polygons
PolExtract <- function(registration) {
  d1 <- do.call(rbind, lapply(1:registration$atlas$numRegions, function(x) {
    d <- data.frame(x = registration$atlas$outlines[[x]]$xlT, y = registration$atlas$outlines[[x]]$ylT)
    d$region <- x
    d$hemisphere <- "left"
    d$region_combined <- paste0(d$region, "_", d$hemisphere)
    d$color <- registration$atlas$col[x] %>% as.character()
    return(d)
  }))
  d2 <- do.call(rbind, lapply(1:registration$atlas$numRegions, function(x) {
    d <- data.frame(x = registration$atlas$outlines[[x]]$xrT, y = registration$atlas$outlines[[x]]$yrT)
    d$region <- x
    d$hemisphere <- "right"
    d$region_combined <- paste0(d$region, "_", d$hemisphere)
    d$color <- registration$atlas$col[x] %>% as.character()
    return(d)
  }))
  df <- rbind(d1, d2)
  return(df)
}

gg <- na.omit(PolExtract(regi))
ggplot() +
  geom_polygon(data = gg, aes(x, y, group = region_combined), fill = gg$color, alpha = 1, color = NA) +
  geom_polygon(data = gg, aes(x, y, group = region_combined), alpha = 0, color = "black") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_reverse(expand = c(0,0)) +
  theme_void()
