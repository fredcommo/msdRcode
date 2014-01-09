# Use results from predict_MSD_CK_2013_10_17 and predict_MSD_RK_2013_10_17
op <- par(no.readonly = TRUE)
require(rGithubClient)
git <- getRepo('fredcommo/proteinArrays')
sourceRepoFile(git, "predict_protArrays_helpers.R")

# Set inDir
inDir <- "path/to/MSDresults"

# Read all the results to merge
ck <- read.csv(file.path(indir, "MSD_Results_CK_2013-10-17.txt"), header = TRUE, sep = "\t")
rk <- read.csv(file.path(indir, "MSD_Results_RK_2013-10-17.txt"), header = TRUE, sep = "\t")

# Filter on samples of interest
ck <- ck[grep("^CK", ck$sample_vs_ref),]
rk <- rk[grep("^RK", rk$sample_vs_ref),]

# Merge
all <- rbind.data.frame(ck, rk)
all <- mergeTables(all)
colnames(all) <- gsub("_(.*)| ", "", colnames(all))
Ids <- gsub("_(.*)", "", rownames(all))
subIds <- gsub("(.*)(R\\d+)(.*)", "\\2", rownames(all))

# Plot
png(file.path(indir, "heatmap_samples_msd_2013-10-17.png"), width = 800, height = 600)
HeatMap(log10(all), Ids, subIds)
dev.off()
