library(plyr)
library(dplyr)
library(readxl)
library(purrr)
library(writexl)


# Paired-pulse recordings should be pre-analyzed in Clampfit. 
# Traces recorded at increasing stimulation intensities should be concatenated and baseline should be normalized. 
# Then, the amplitude of first (P1) and second peak (P2) should be calculated.


# This script will help analyze the recordings for each cell, computing the 
# estimated number of CF inputs to a given PC and their relative strength. 
# The following files will be generated:
  # - analysis_DI_DR.RDS and analysis_DI_DR.xlsx: dataframes with the number of estimated CF inputs,
  #                                               their amplitude (pA), DI, DR, PPR per input and average PPR.
  # - ls_subset.RDS: list of dfs containing the sweeps assigned to each input.
  # - cell 24_densityplot.png: summary plot showing the density distribution and where the amplitudes were subset.


# 1) Crate analysis directory ####
cell = "cell 4"
analysis_path <- paste0(cell, "_CF inputs analysis")
dir.create(analysis_path, showWarnings = FALSE)


# 2) Load detected amplitudes from Clampfit####
# All the calculated values are in a single .xlsx file, with one sheet per cell
amp <- read_xlsx("Peak detection.xlsx", sheet = cell)  %>% as.data.frame(.)



# 3) Remove facilitated pulses (P2 > P1) ####
# Calculate PPR
amp$ppr <- amp$p2/amp$p1

# Remove sweeps with facilitated pulses or failed stimulation
amp <- amp[!((amp$ppr>1)| (amp$ppr < 0)),] 

# Remove positive values traces (normalization artifacts)
amp <- amp[!((amp$p1 > 0) | (amp$p2 > 0)),]



# 4) Assign peaks to CF inputs ####
# Compute density distribution of first peak and calculate local minima
p1 <- amp$p1
d <- density(p1, kernel = "gaussian", bw = "SJ")

localmin_idx <- which(diff(sign(diff(d$y))) > 0) + 1L
localmin <- d$x[localmin_idx]

# Assign each value to an input based on density distirbution minima
breaks <- c(-Inf, localmin, Inf)
amp$interval_id <- cut(amp$p1, breaks = breaks, labels = FALSE, include.lowest = TRUE, right = TRUE)

# Assign each point of the density curve to the same intervals
den <- data.frame(dx = d$x, dy = d$y)
den$interval_id <- cut(den$dx, breaks = breaks, labels = FALSE, include.lowest = TRUE, right = TRUE)

# Summarize each density / amplitude interval
interval_ids <- seq_len(length(breaks) - 1L)
thresden <- sd(d$y) / 5

interval_summary <- data.frame(interval_id = interval_ids, left = c(min(d$x), localmin), right = c(localmin, max(d$x)), 
                               n_sweeps = tabulate(amp$interval_id,nbins = length(interval_ids)),
                               max_density = vapply(interval_ids, function(i) max(den$dy[den$interval_id == i],
                                na.rm = TRUE), numeric(1)))



## Refine intervals - reject inputs with few sweeps ####
# Choose method to reject inputs:
# filter_method can be "n_sweeps", to remove inuts with 2 seeps or fewer,
# or "density", to exclude inputs whose density does not reach the defined threshold.
filter_method <- "n_sweeps"

# Define minimum number of sweeps to be consider as an independent input
min_sweeps <- 3L                  

valid_interval <- switch(filter_method, n_sweeps = interval_summary$n_sweeps >= min_sweeps,
                         density = interval_summary$max_density >= thresden,
                         stop("filter_method must be 'n_sweeps' or 'density'"))


# Subset with corrected intervals
ls_subset_all <- split(amp, amp$interval_id)
valid_ids <- as.character(interval_summary$interval_id[valid_interval])
ls_subset_inter2 <- ls_subset_all[intersect(valid_ids, names(ls_subset_all))]

ls_subset_inter2 <- lapply(ls_subset_inter2, function(x) {x$interval_id <- NULL
                                                          x})
names(ls_subset_inter2) <- NULL

kept_ranges <- interval_summary[valid_interval, c("left", "right"), drop = FALSE]
intervals2 <- sort(unique(c(kept_ranges$left, kept_ranges$right)))


## Plot and save ####
densityfile <- paste0(analysis_path, "/", cell, "_densityplot.png")
png(file = densityfile, width = 2000, height = 800)
plot(d)
abline(h = thresden, col = "red")
rect(ybottom = 0, ytop = thresden, xleft = min(d$x), xright = max(d$x), border = NA, col = adjustcolor("red", alpha.f = 0.2))
abline(v = localmin, col = "red")
abline(v = intervals2, col = "green")
dev.off()

saveRDS(ls_subset_inter2, paste0(analysis_path, "/ls_subset.RDS"))
write_xlsx(interval_summary, paste0(analysis_path, "/interval_summary.xlsx"))



# 5) Average sweeps assigned to the same CF input ####
averaged_subsets <- lapply(ls_subset_inter2, colMeans) %>% lapply(., as.data.frame)
averaged_subsets_df <- list_cbind(averaged_subsets)
averaged_subsets_df <- as.data.frame(t(averaged_subsets_df))
averaged_subsets_df <- averaged_subsets_df[,c(3,5,6)]



# 6) Analyze relative strenght between inputs (DI, DR) ####
## Calculate real inputs amplitude ####
# Since amplitude steps are the result of the progressive recruitment of more CFs
# when increasing stimulation intensity, bigger amplitudes in multi-innervated
# cells are the result of the combined amplitude of several inputs. Therefore,
# it is necessary to subtract it to access the real input value
averaged_subsets_df <- averaged_subsets_df[order(averaged_subsets_df$p1),]
p1 <- averaged_subsets_df[,1]
p2 <- averaged_subsets_df[,2]
ncf <- length(p1)

p1_subs <- c()
p2_subs <- c()
for (val in 1:ncf) {
  p1_subs[val] <- p1[val] - p1[val+1]
  p2_subs[val] <- p2[val] - p2[val+1]
}

p1_subs[ncf] <- p1[ncf] 
p2_subs[ncf] <- p2[ncf]

# Bind new substracted inputs into a dataframe
avinputs_subs <- cbind(p1_subs, p2_subs) %>% as.data.frame(.) 
colnames(avinputs_subs) <- c("P1", "P2")
avinputs_subs = avinputs_subs[order(-avinputs_subs$P1),]  


## Calculate PPR per averaged input and per cell ####
avinputs_subs$PPRsub <- avinputs_subs$P2 / avinputs_subs$P1
avinputs_subs$PPRav <- mean(avinputs_subs$PPRsub)


## Calculate disparity index (DI = -S.D./mean) ####
m = mean(avinputs_subs$P1)
sd = sd(avinputs_subs$P1)
DI = -(sd(avinputs_subs$P1))/(mean(avinputs_subs$P1)) #Calculates DI
avinputs_subs$DI <- DI


## Calculate disparity ratio (DR = (E(Ai/An))/(N-1) ) ####
ai <- avinputs_subs$P1 %>% as.numeric(.)

ratio <- as.numeric()
n <- length(ai)
for (a in 1:length(ai)) { 
  ratio[a] <- ((ai[a])/(ai[n])) #Devides each value by the biggest amplitude
}

ratio = ratio[-n] 
DR = sum(ratio) / (n-1) 
avinputs_subs$DR <- DR


## Add number of CF inputs ####
avinputs_subs$N_CF <- nrow(avinputs_subs)



# Save results ####
saveRDS(avinputs_subs, paste0(analysis_path, "/analysis_DI_DR.RDS", sep=""))
write_xlsx(avinputs_subs, path = paste0(analysis_path, "/analysis_DI_DR.xlsx", sep=""))

