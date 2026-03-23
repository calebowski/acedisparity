# Test script to verify environment and paths
cat("Starting R script...\n")

# Try loading the core libraries
libs <- c("dplyr", "tidyr", "lme4", "lmerTest")
for (lib in libs) {
  if (require(lib, character.only = TRUE)) {
    cat(paste("Successfully loaded:", lib, "\n"))
  } else {
    cat(paste("FAILED to load:", lib, "\n"))
  }
}

# Test folder access
test_path <- "/mnt/parscratch/users/bip24cns/acedisparity/continuous/"
if (dir.exists(test_path)) {
  cat(paste("Success: Access to parscratch confirmed at", test_path, "\n"))
} else {
  cat("ERROR: Cannot access parscratch path\n")
}

cat("Test script complete.\n")