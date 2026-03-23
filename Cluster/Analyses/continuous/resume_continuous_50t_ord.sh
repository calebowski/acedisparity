#!/bin/bash
#SBATCH --job-name=RESUME_acecontinuous50tord
#SBATCH --output=/users/bip24cns/acedisparity/continuous/logs/resume_ord_continuous_50t_%A_%a.out
#SBATCH --error=/users/bip24cns/acedisparity/continuous/logs/resume_ord_continuous_50t_%A_%a.err
#SBATCH --array=18,26,52,60,90,99       
#SBATCH --nodes=1                     # 1 node per task
#SBATCH --ntasks=5                    # 1 R process per task
#SBATCH --cpus-per-task=5            # adjust if using parallel inside R
#SBATCH --mem=16G                    # adjust based on memory needs
#SBATCH --time=65:00:00               
#SBATCH --mail-user=cnscutt1@sheffield.ac.uk
#SBATCH --mail-type=ALL

module load R/4.4.1-foss-2022b

export R_LIBS_USER=/users/$USER/R/x86_64-pc-linux-gnu-library/4.4


REPLICATE=${SLURM_ARRAY_TASK_ID}
MODELS=("bm" "bm_t" "ou_w" "ou_st" "ou_sh")

Rscript /users/$USER/acedisparity/continuous/scripts/ace_continuous_50t_generate_tree.R ${REPLICATE}

for i in "${!MODELS[@]}"; do
    MODEL="${MODELS[$i]}"
    srun --ntasks=1 --cpus-per-task=5 --exclusive --export=ALL \
        Rscript /users/$USER/acedisparity/continuous/scripts/ace_continuous_process_model.R ${REPLICATE} ${MODEL} &
done

wait

Rscript /users/$USER/acedisparity/continuous/scripts/combine_save_continuous_models.R ${REPLICATE} "50t"

echo "Completed replicate ${REPLICATE}"