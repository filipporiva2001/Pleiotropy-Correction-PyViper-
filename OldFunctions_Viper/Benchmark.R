setwd("/Users/friva/Library/CloudStorage/OneDrive-SanfordBurnhamPrebysMedicalDiscoveryInstitute/Desktop/SBP_CODE/PleiotropyProject/BRCA")

## ---- Config: set your file paths ----
r_csv <- "NES_RVIPER_precorrectionBRCA1.csv"
p_csv <- "NES_pyVIPER_precorrection_BRCA_github.csv"

## ---- Helper: read a matrix allowing first column as rownames ----
read_named_matrix <- function(path) {
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  # If first column looks like row identifiers (non-numeric), use it as rownames
  if (!is.numeric(df[[1]]) && !anyDuplicated(df[[1]])) {
    rn <- df[[1]]
    df <- df[ , -1, drop = FALSE]
    rownames(df) <- rn
  }
  # Convert to numeric matrix
  m <- as.matrix(df)
  storage.mode(m) <- "numeric"
  return(m)
}

R_mat <- read_named_matrix(r_csv)  # R-based VIPER NES
P_mat <- read_named_matrix(p_csv)  # Python pyVIPER NES

## ---- Align matrices by intersecting row/column names ----
common_rows <- intersect(rownames(R_mat), rownames(P_mat))
common_cols <- intersect(colnames(R_mat), colnames(P_mat))

if (length(common_rows) == 0 || length(common_cols) == 0) {
  stop("No overlapping rows or columns between the two matrices after name matching.")
}

R_aln <- R_mat[common_rows, common_cols, drop = FALSE]
P_aln <- P_mat[common_rows, common_cols, drop = FALSE]

## ---- Sanity checks ----
stopifnot(identical(rownames(R_aln), rownames(P_aln)))
stopifnot(identical(colnames(R_aln), colnames(P_aln)))

## ---- Compute differences ----
diff_abs <- abs(R_aln - P_aln)           # |R_ij - P_ij|
denom    <- abs(R_aln)                   # |R_ij|

## Exclude positions where denom == 0 to avoid division by zero for MRD
mask <- denom > 0 & is.finite(denom)
num_excluded <- sum(!mask)
total_entries <- length(denom)

## Mean Relative Difference (%) over valid entries
mrd_percent <- mean((diff_abs[mask] / denom[mask]) * 100, na.rm = TRUE)

## Max Absolute Difference over all aligned entries
max_abs_diff <- max(diff_abs, na.rm = TRUE)

## ---- Report ----
cat(sprintf("Aligned size: %d rows x %d cols\n",
            nrow(R_aln), ncol(R_aln)))
cat(sprintf("Entries excluded from MRD due to |R_ij| = 0: %d of %d (%.2f%%)\n",
            num_excluded, total_entries, 100 * num_excluded / total_entries))
cat(sprintf("Mean Relative Difference (%%): %.6f\n", mrd_percent))
cat(sprintf("Max Absolute Difference: %.6f\n", max_abs_diff))
## ---- Find TF (row) and sample (column) with max absolute difference ----
max_pos <- which(diff_abs == max_abs_diff, arr.ind = TRUE)

# Extract the corresponding TF name and sample name
max_tf <- rownames(diff_abs)[max_pos[1, "row"]]
max_sample <- colnames(diff_abs)[max_pos[1, "col"]]

# Also report their R and Python NES values
R_value <- R_aln[max_tf, max_sample]
P_value <- P_aln[max_tf, max_sample]

cat("\n--- Maximum Absolute Difference Details ---\n")
cat(sprintf("TF with max abs difference: %s\n", max_tf))
cat(sprintf("Sample / Column: %s\n", max_sample))
cat(sprintf("R-VIPER NES value: %.6f\n", R_value))
cat(sprintf("pyVIPER NES value: %.6f\n", P_value))
cat(sprintf("Absolute Difference: %.6f\n", max_abs_diff))
