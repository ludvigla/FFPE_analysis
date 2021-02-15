library(ggplot2)
library(patchwork)

dirs <- paste0("~/FFPE/", c("200820_lung_human_organoid", "200924_kidney_human_organoid", 
          "201016_brain_mouse_NA", "201016_lung_human_Covid19",
          "201026_cancer_human_archival"))

trim.files <- unlist(sapply(dirs, function(d) {
  list.files(path = d, pattern = "cutadapt_log.txt", recursive = TRUE, full.names = TRUE)
}))

samplenames <- do.call(rbind, strsplit(trim.files, "/"))[, 6]
dates <- do.call(rbind, strsplit(trim.files, "/"))[, 5]

read.tso.stats <- function(f) {
  TSO <- read.table(f, header = T, skip = 22, nrows = 29, sep = "\t")
  return(TSO)
}
read.polyA.stats <- function(f) {
  polyA <- read.table(f, header = T, skip = 69, nrows = length(readLines(f)) - 69, sep = "\t")
  return(polyA)
}
read.tot.count <- function(f) {
  L <- readLines(f)[8]
  count <- as.numeric(paste(unlist(regmatches(L, gregexpr("[[:digit:]]+", L))), collapse = ""))
  return(count)
}

# summarize TSO trimming stats
tso.df <- do.call(rbind, lapply(seq_along(trim.files), function(i) {
  f <- trim.files[i]
  TSO <- read.tso.stats(f)
  tot.count <- read.tot.count(f)
  TSO$fraction <- TSO$count/tot.count
  TSO$sample <- samplenames[i]
  TSO$date <- dates[i]
  return(TSO)
}))

p1 <- ggplot(tso.df, aes(length, fraction, color = sample)) +
  geom_line() +
  geom_point(stat = "identity", size = 1, color = "black") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#114477", "#4477AA", "#77AADD", "#117755", "#44AA88", "#99CCBB", "#777711", "#AAAA44", "#DDDD77", "#771111", "#AA4444", "#DD7777", "#771144", "#AA4477", "#DD77AA")) +
  facet_grid(~date) +
  guides(color = F) +
  labs(x = "Adapter length", y = "fraction of total reads") +
  theme_minimal()


# summarize polyA trimming stats
polyA.df <- do.call(rbind, lapply(seq_along(trim.files), function(i) {
  f <- trim.files[i]
  polyA <- read.polyA.stats(f)
  tot.count <- read.tot.count(f)
  polyA$fraction <- polyA$count/tot.count
  polyA$sample <- samplenames[i]
  polyA$date <- dates[i]
  return(polyA)
}))

p2 <- ggplot(polyA.df, aes(length, fraction, color = sample)) +
  geom_line() +
  geom_point(stat = "identity", size = 1, color = "black") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#114477", "#4477AA", "#77AADD", "#117755", "#44AA88", "#99CCBB", "#777711", "#AAAA44", "#DDDD77", "#771111", "#AA4444", "#DD7777", "#771144", "#AA4477", "#DD77AA")) +
  facet_grid(~date) +
  guides(color = F) +
  labs(x = "polyA length", y = "fraction of total reads") +
  theme_minimal()

jpeg(filename = "~/FFPE/analysis/trimming_stats.jpeg", width = 3500, height = 2500, res = 300)
print(p1 / p2)
dev.off()