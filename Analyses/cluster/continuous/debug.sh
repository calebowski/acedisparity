#!/bin/bash
# Save this as debug_tree_issue.sh and run: bash debug_tree_issue.sh

echo "=== FILESYSTEM DEBUGGING SCRIPT ==="
echo "Date: $(date)"
echo "User: $(whoami)"
echo "Node: $(hostname)"
echo

# Set variables to match your script
BASE_PATH="/mnt/parscratch/users/bip24cns/acedisparity/continuous/100t"
TREES_DIR="${BASE_PATH}/trees"
JOB_ID="8453504"
REPLICATE_ID=24

echo "=== TESTING FILESYSTEM ACCESS ==="
echo "Base path: $BASE_PATH"
echo "Trees dir: $TREES_DIR"
echo

# Test base path
echo "1. Base path tests:"
echo "   Exists: $(test -d "$BASE_PATH" && echo "YES" || echo "NO")"
echo "   Readable: $(test -r "$BASE_PATH" && echo "YES" || echo "NO")"
echo "   Writable: $(test -w "$BASE_PATH" && echo "YES" || echo "NO")"
echo "   Executable: $(test -x "$BASE_PATH" && echo "YES" || echo "NO")"
echo

# Test parent directory permissions
PARENT_DIR=$(dirname "$TREES_DIR")
echo "2. Parent directory tests ($PARENT_DIR):"
echo "   Exists: $(test -d "$PARENT_DIR" && echo "YES" || echo "NO")"
echo "   Readable: $(test -r "$PARENT_DIR" && echo "YES" || echo "NO")"
echo "   Writable: $(test -w "$PARENT_DIR" && echo "YES" || echo "NO")"
echo "   Executable: $(test -x "$PARENT_DIR" && echo "YES" || echo "NO")"
echo

# Check if trees directory exists
echo "3. Trees directory tests:"
echo "   Exists: $(test -d "$TREES_DIR" && echo "YES" || echo "NO")"
if [ -d "$TREES_DIR" ]; then
    echo "   Readable: $(test -r "$TREES_DIR" && echo "YES" || echo "NO")"
    echo "   Writable: $(test -w "$TREES_DIR" && echo "YES" || echo "NO")"
    echo "   Executable: $(test -x "$TREES_DIR" && echo "YES" || echo "NO")"
fi
echo

# Test creating trees directory
echo "4. Testing directory creation:"
if [ ! -d "$TREES_DIR" ]; then
    echo "   Attempting to create $TREES_DIR..."
    if mkdir -p "$TREES_DIR" 2>/dev/null; then
        echo "   SUCCESS: Directory created"
        echo "   Now writable: $(test -w "$TREES_DIR" && echo "YES" || echo "NO")"
    else
        echo "   FAILED: Could not create directory"
        echo "   Error details:"
        mkdir -p "$TREES_DIR"
    fi
else
    echo "   Directory already exists"
fi
echo

# Test file creation
echo "5. Testing file creation:"
TEST_FILE="${TREES_DIR}/test_file_${REPLICATE_ID}.txt"
echo "   Attempting to create test file: $TEST_FILE"
if echo "test content" > "$TEST_FILE" 2>/dev/null; then
    echo "   SUCCESS: File created"
    echo "   File size: $(stat -c%s "$TEST_FILE" 2>/dev/null || echo "unknown") bytes"
    rm -f "$TEST_FILE" 2>/dev/null
    echo "   Cleanup: $(test -f "$TEST_FILE" && echo "FAILED" || echo "SUCCESS")"
else
    echo "   FAILED: Could not create file"
    echo "   Error details:"
    echo "test content" > "$TEST_FILE"
fi
echo

# Check disk space and quotas
echo "6. Disk space and quotas:"
echo "   Disk usage for $BASE_PATH:"
df -h "$BASE_PATH"
echo
echo "   User quotas (if available):"
quota -u $(whoami) 2>/dev/null || echo "   No quota information available"
echo

# Check filesystem type and mount options
echo "7. Filesystem information:"
echo "   Mount point for $BASE_PATH:"
df "$BASE_PATH" | tail -n 1
echo
echo "   Filesystem type and options:"
mount | grep "$(df "$BASE_PATH" | tail -n 1 | awk '{print $1}')" || echo "   Could not determine mount options"
echo

# List existing files
echo "8. Existing files in base directory:"
ls -la "$BASE_PATH" 2>/dev/null || echo "   Could not list base directory"
echo
if [ -d "$TREES_DIR" ]; then
    echo "   Files in trees directory:"
    ls -la "$TREES_DIR" 2>/dev/null || echo "   Could not list trees directory"
fi
echo

# Check for any existing tree files
echo "9. Checking for existing tree files:"
find "$BASE_PATH" -name "*tree*" -type f 2>/dev/null | head -10 || echo "   No tree files found or permission denied"
echo

echo "=== R DEBUGGING SCRIPT ==="
echo "Now testing with R..."

# Create R debugging script
cat > debug_tree.R << 'EOF'
args <- c("24")  # Hardcode replicate 24 for testing
replicate_id <- as.numeric(args[1])

cat("=== R DEBUGGING SESSION ===\n")
cat("R version:", R.version.string, "\n")
cat("Working directory:", getwd(), "\n")
cat("User:", Sys.getenv("USER"), "\n")
cat("\n")

# Test basic R file operations
base_path <- "/mnt/parscratch/users/bip24cns/acedisparity/continuous/100t/"
job_id <- 8453504

write.path <- function(subfolder, filename) {
  paste0(base_path, subfolder, "/", job_id, "_", sprintf(filename, replicate_id))
}

# Test path construction
tree_path <- write.path("trees", "tree_%03d.tre")
trees_dir <- dirname(tree_path)

cat("Constructed paths:\n")
cat("  Base path:", base_path, "\n")
cat("  Tree path:", tree_path, "\n")
cat("  Trees dir:", trees_dir, "\n")
cat("\n")

# Test R file operations
cat("R file operation tests:\n")
cat("  Base path exists:", dir.exists(base_path), "\n")
cat("  Trees dir exists:", dir.exists(trees_dir), "\n")
cat("  Base path access:", file.access(base_path, mode = 2), "\n")
cat("  Trees dir access:", if(dir.exists(trees_dir)) file.access(trees_dir, mode = 2) else "N/A", "\n")
cat("\n")

# Test directory creation
cat("Testing directory creation with R:\n")
if(!dir.exists(trees_dir)) {
  cat("  Attempting to create:", trees_dir, "\n")
  result <- tryCatch({
    dir.create(trees_dir, recursive = TRUE, showWarnings = TRUE)
  }, error = function(e) {
    cat("  ERROR:", e$message, "\n")
    return(FALSE)
  })
  cat("  Creation result:", result, "\n")
  cat("  Directory exists after creation:", dir.exists(trees_dir), "\n")
}

# Test file writing
cat("Testing file writing:\n")
test_content <- "This is a test file"
test_file <- file.path(trees_dir, "test_r_file.txt")

write_result <- tryCatch({
  writeLines(test_content, test_file)
  cat("  SUCCESS: Test file written\n")
  cat("  File exists:", file.exists(test_file), "\n")
  cat("  File size:", file.size(test_file), "bytes\n")
  
  # Cleanup
  unlink(test_file)
  cat("  Cleanup successful:", !file.exists(test_file), "\n")
  TRUE
}, error = function(e) {
  cat("  ERROR writing test file:", e$message, "\n")
  FALSE
})

cat("\n=== END R DEBUGGING ===\n")
EOF

# Run the R debugging script
echo "Running R debugging script..."
module load R/4.4.1-foss-2022b 2>/dev/null || echo "Could not load R module"
Rscript debug_tree.R

echo
echo "=== DEBUGGING COMPLETE ==="
echo "Check the output above for any permission or filesystem issues."