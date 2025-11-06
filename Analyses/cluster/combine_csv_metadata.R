library(data.table)
meta_dir <- "/mnt/parscratch/users/bip24cns/acedisparity/trees/overallDisparity/metadata"
files <- list.files(meta_dir, pattern="^metadata_\\d{3}\\.csv$", full.names=TRUE)
all_meta <- data.table::rbindlist(lapply(sort(files), data.table::fread), fill = TRUE)
out_file <- file.path(dirname(meta_dir), "metadata_all.csv")
tmp <- tempfile(pattern="meta_all_", tmpdir=dirname(meta_dir))
data.table::fwrite(all_meta, tmp)
file.rename(tmp, out_file)
cat("Wrote combined metadata to", out_file, "\n")