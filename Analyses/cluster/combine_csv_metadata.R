meta_dir <- "/home/caleb/Documents/PhD/acedisparity/Data/cluster/trees/large/metadata/"
files <- list.files(meta_dir, pattern = "^metadata_\\d{3}\\.csv$", full.names = TRUE)

# Read all CSV files and combine
all_meta <- do.call(rbind, lapply(sort(files), function(f) {
  read.csv(f, stringsAsFactors = FALSE)
}))

# Write to temp file then rename (atomic write)
out_file <- file.path(dirname(meta_dir), "metadata_all.csv")
tmp <- tempfile(pattern = "meta_all_", tmpdir = dirname(meta_dir), fileext = ".csv")
write.csv(all_meta, tmp, row.names = FALSE, quote = TRUE)
file.rename(tmp, out_file)
cat("Wrote combined metadata to", out_file, "\n")