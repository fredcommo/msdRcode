# Predict - MSD data Rcode
op <- par(no.readonly = TRUE)
require(preprocessCore)
require(rGithubClient)
git <- getRepo('fredcommo/msdRcode')
sourceRepoFile(git, "msd_helper_functions.R")

op <- par(no.readonly = TRUE)
dev.off()
today <- format(Sys.Date(), format="%Y_%m_%d")
subTit <- paste('Meso Scale discovery -', today)
  
# Set your inDir (containing data)
#inDir <- "/Users/fredcommo/Documents/Predict/Kinome_Swanton_MSD"
inDir <- "path/to/MSD_data"
outDir <- file.path(inDir, sprintf("Results_%s", today))
dir.create(outDir)

# Choose the file you want to work with
fileName <- "All Data_MSD_template_final_RK_2013_10_14.txt"
outputName <- sprintf("%s_Results_%s", gsub(".txt", "", fileName), today)
rawData <- read.csv(file.path(inDir, fileName), sep = '\t')
head(rawData)

# At a glance...
boxplot(log(rawData$value) ~ rawData$sample, par(mar = c(10, 5, 4, 2),las = 3))

par(mfrow=c(1, 2), mar = c(10, 5, 4, 2), cex.lab=1.5, las = 3)
  # Human serum profiles
humSer <- rawData[grep("human_serum", rawData$sample),]
boxplot(humSer$value ~ humSer$target, xlab="", ylab = "Signal", main = "RK24: Human serum")
title(xlab = "proteins", mgp=c(7.5, 1, 0))
  # GAPDH profiles
gapdh <- rawData[grep("^GAPDH", rawData$target),]
boxplot(gapdh$value ~ gapdh$sample, xlab = "", main="RK24: GAPDH")
title(xlab = "samples", mgp=c(7.5, 1, 0))
par(op)

########################################

# Adjust and correct signal
adjData <- adjustHS(rawData)
adjData <- adjustGAPDH(adjData)

# Negative values are replaced by random positive values
if(any(adjData$value <= 0)){
  idx <- which(adjData$value <= 0)
  adjData$value[idx] <- 1 + abs(rnorm(length(idx), 0, 1))
}

# Run ratios on samples: all the samples may not have all the same proteins.
A431ratios <- rbind.data.frame(.ratio(adjData, "A431_Ctrl", "A431_EGF", "mTOR", "p-mTOR"),
                               .ratio(adjData, "A431_Ctrl", "A431_EGF", "Akt", "p-Akt (ser473)"),
                               .ratio(adjData, "A431_Ctrl", "A431_EGF", "Akt", "p-Akt (thr308)"))

MCFratios <- .ratio(adjData, "MCF-7_Ctrl", "MCF-7_Insulin", "p70S6K", "p-p70S6K")

CKratios <- lapply(c("RK24_R1", "RK24_R3", "RK24_R4", "RK24_R9"), function(samp){
  cat(samp, '\n')
  ref <- "RK24_N"
  ratios <- rbind.data.frame(.ratio(adjData, ref, samp, "mTOR", "p-mTOR"),
                             .ratio(adjData, ref, samp, "Akt", "p-Akt (ser473)"),
                             .ratio(adjData, ref, samp, "Akt", "p-Akt (thr308)"),
                             .ratio(adjData, ref, samp, "p70S6K", "p-p70S6K"),
                             .ratio(adjData, ref, samp, "4E-BP1", "p-4E-BP1"))
  return(ratios)
})
CKratios <- as.data.frame(do.call(rbind, CKratios))

# Merge Results
Results <- rbind.data.frame(A431ratios, MCFratios, CKratios)

# Plot Results
png(file = file.path(outDir, sprintf("%s_%s.png", outputName, today)), width = 1000, height = 800)
par(las = 3, mar = c(7, 6, 4, 2))
n <- nlevels(Results$target)
Cols <- grey(seq(1/n, 1, len = n))
barplot(log2(Results$RR), name = gsub("_vs(.*)", "", Results$sample_vs_ref),
        col = Cols[Results$target], ylab = expression(log[2](RR)), main = subTit)
legend("topright", legend = unique(Results$target), fill = Cols)
par(op)
dev.off()

# Save results
write.table(Results, file.path(outDir, sprintf("%s_%s.txt", outputName, today)), sep = '\t', row.names = FALSE)


# Old code
# for(i in c(2, 4)){
#   png(file = paste0(getwd(), '/MSD_', colnames(ratios)[i+1], '_', today, '.png'), width = 1600, height = 1000)
#   par(mfrow = c(1, 2), mar = c(15, 7, 4, 2), cex.axis = 1.5, cex.lab = 1.75, cex.main = 2.5)
#   barplot(ratios[,i], main = colnames(ratios)[i], ylab = 'ratio(Total/Phospo)',
#           ylim = range(ratios[,i:(i+1)], na.rm = TRUE),
#           names.arg = ratios$prot, las = 3)
#   barplot(ratios[,(i+1)], main = colnames(ratios)[i+1],
#           ylim = range(ratios[,i:(i+1)], na.rm = TRUE),
#           names.arg = ratios$prot, las = 3)
#   par(op)
#   dev.off()
# }

# write.table(ratios, 'MSD_PA_86.xls', sep = '\t')

