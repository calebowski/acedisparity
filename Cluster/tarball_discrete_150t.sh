#!/bin/bash
#SBATCH --job-name=pack_folders_discrete_150t
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=96:00:00
#SBATCH --mem=10G
#SBATCH --mail-user=cnscutt1@sheffield.ac.uk
#SBATCH --mail-type=ALL

echo "Starting the automated folder packaging job..."

# 1. Navigate to the parent directory containing all your data folders
cd /mnt/parscratch/users/bip24cns/acedisparity/discrete/crown/150t/

# tar -cf 50t_continuous.tar 50t

# 2. Loop through every item in this directory that is a folder (indicated by */)
for folder in */; do
    
    # 3. Strip the trailing slash off the folder name 
    # (e.g., turns '50t/' into a clean '50t' so our file isn't named '50t/.tar')
    clean_name="${folder%/}"
    
    echo "Starting to pack: $clean_name -> discrete_150t_${clean_name}.tar"
    
    # 4. Run the tar command (using -cf for Create File, dropping the verbose output)
    tar -cf "discrete_150t_${clean_name}.tar" "$clean_name"
    
    echo "Finished packing: $clean_name"
    echo "---------------------------------------------------"
    
done

echo "All folders have been successfully packaged!"



