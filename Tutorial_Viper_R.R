# ------------------------------------------------------------
## Libraries
## ------------------------------------------------------------
library(viper)
library(Biobase)

## ------------------------------------------------------------
# Set up directory
## ------------------------------------------------------------
setwd("/Users/friva/Library/CloudStorage/OneDrive-SanfordBurnhamPrebysMedicalDiscoveryInstitute/Desktop/SBP_CODE/PleiotropyProject/BRCA")

## ------------------------------------------------------------
## Helper: table -> regulon (same shape viper expects)
## ------------------------------------------------------------
TableToInteractome <- function(net_table){
  blank_reg <- list(c(), c()); names(blank_reg) <- c("tfmode", "likelihood")
  my_reg <- list()
  pb <- txtProgressBar(min = 0, max = nrow(net_table), initial = 0, style = 3)
  for(i in 1:nrow(net_table)){
    my_reg[[net_table[i,"regulator"]]] <- blank_reg
    setTxtProgressBar(pb, i)
  }
  close(pb)
  pb <- txtProgressBar(min = 0, max = nrow(net_table), initial = 0, style = 3)
  for(i in 1:nrow(net_table)){
    new_w <- net_table[i,"likelihood"]; new_m <- net_table[i,"mor"]
    names(new_w) <- net_table[i,"target"]; names(new_m) <- net_table[i,"target"]
    my_reg[[net_table[i,"regulator"]]][["likelihood"]] <- c(my_reg[[net_table[i,"regulator"]]][["likelihood"]], new_w)
    my_reg[[net_table[i,"regulator"]]][["tfmode"]]     <- c(my_reg[[net_table[i,"regulator"]]][["tfmode"]],     new_m)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  my_reg
}

## ------------------------------------------------------------
## 1) Load data
## ------------------------------------------------------------
expr <- read.delim("DEG-BRCA.tsv", row.names = 1, check.names=FALSE)
expr_matrix <- t(as.matrix(expr))  # samples in columns
dset <- ExpressionSet(assayData = expr_matrix)

interactome <- read.delim("Interactome_prunedBRCA1.tsv", stringsAsFactors = FALSE)
regulon <- TableToInteractome(interactome)


## ------------------------------------------------------------
## 2) VIPER pre-pleiotropy NES
## ------------------------------------------------------------
tf_nes_original <- viper(
  exprs(dset), regulon,
  pleiotropy = FALSE,
  method = "none",
  nes = TRUE,
  eset.filter = TRUE,
  minsize = 0,    #disables regulon-size filtering
  cores = 1,
  verbose = TRUE
)
write.csv(tf_nes_original, "NES_RVIPER_precorrectionBRCA1.csv", row.names = TRUE)

dd_default <- viperSimilarity(tf_nes_original)              # returns a signatureDistance (square matrix)

dd_default[1:5, 1:5]
scale(dd_default)[1:5, 1:5]
as.matrix(as.dist(dd_default))[1:5, 1:5]


getAnywhere(viperSimilarity)
