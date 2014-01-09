# MSD plates analysis
# Helper functions
# FC 2014_01_08

require(rGithubClient)
git <- getRepo('fredcommo/proteinArrays')
sourceRepoFile(git, "predict_protArrays_helpers.R")

# compute Ratios
adjustGAPDH <- function(Table){
  gapdh <- Table[which(Table$target == "GAPDH"),]
  gapdh <- by(gapdh$value, gapdh$sample, mean, na.rm = TRUE)
  alpha <- as.vector(max(gapdh, na.rm = TRUE)/gapdh)
  names(alpha) <- names(gapdh)
  for(id in names(alpha)){
    a <- alpha[which(names(alpha) == id)]
    if(!is.na(a))
      Table$value[which(Table$sample == id)] <- Table$value[which(Table$sample == id)]*a
  }
  return(Table)
}
adjustHS <- function(Table){
  targets <- levels(Table$target)
  out <- sapply(targets, function(target){
    x <- Table$value[which(Table$target == target)]
    h <- Table$value[which(Table$target == target & Table$sample == "human_serum")]
    if(length(h) == 0)
      return(x)
    else
      return(x - mean(h, na.rm = TRUE))
  })
  for(target in names(out))
    Table$value[which(Table$target == target)] <- out[[which(names(out) == target)]]
  return(Table)
}
.ratio <- function(Table, sample1, sample2, target, pTarget){
  x <- Table$value[which(Table$target == target & Table$sample == sample1)]
  px <- Table$value[which(Table$target == pTarget & Table$sample == sample1)]
  y <- Table$value[which(Table$target == target & Table$sample == sample2)]
  py <- Table$value[which(Table$target == pTarget & Table$sample == sample2)]
  samp1 <- as.vector(sapply(1:length(x), function(i, j){px[j]/x[i]}, 1:length(px)))
  samp2 <- as.vector(sapply(1:length(y), function(i, j){py[j]/y[i]}, 1:length(py)))
  test <- Resamp.T.v4(log2(samp1), log2(samp2))
  Mx <- 2^test$Meanx
  My <- 2^test$Meany
  output <- data.frame("sample_vs_ref" = sprintf("%s_vs_%s", sample2, sample1),
                       "target" = sprintf("%s_vs_%s", pTarget, target),
                       My, Mx, RR = My/Mx, "p" = test$pValue)
  colnames(output)[3:4] <- c("sample", "reference")
  return(as.data.frame(output))
}
mergeTables <- function(msd){
  split.msd <- split(msd, as.character(msd$sample_vs_ref))
  out <- lapply(split.msd, function(x)x$RR)
  out <- as.data.frame(do.call(rbind, out))
  colnames(out) <- split.msd[[1]]$target
  return(out)
}
