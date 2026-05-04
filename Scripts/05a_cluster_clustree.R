# 05a_cluster_clustree.R
# Requires: sample_name  (from 00_setup.R)
# ─────────────────────────────────────────────────────────────────

# Load the Pre-processed dataset
obj <- readRDS(paste0(sample_name, "_post_banksy.rds"))
cat("Loaded:", paste0(sample_name, "_post_banksy.rds"), "\n")

# Prompt: enter 5 resolutions
repeat {
  cat("\nEnter 5 cluster resolutions between 0 and 1.5.\n")
  res_input   <- trimws(readline(prompt = "Comma-separated (e.g. 0.1,0.3,0.5,0.8,1.2): "))
  resolutions <- suppressWarnings(as.numeric(strsplit(res_input, ",")[[1]]))
  if (length(resolutions) == 5 &&
      !any(is.na(resolutions))  &&
      all(resolutions >= 0)     &&
      all(resolutions <= 1.5)) break
  cat("Invalid — please enter exactly 5 numbers between 0 and 1.5, comma-separated.\n")
}

cat("\nRunning clustering at:",
    paste(round(resolutions, 2), collapse = ", "), "\n\n")