#!/bin/bash
#SBATCH --job-name=dna_decode_experiments
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=06:00:00
#SBATCH --mem=128G
#SBATCH --output=experiments_result_%j.log
#SBATCH --error=experiments_error_%j.log

# TSUBAMEでの動的枝刈り+OpenMP並列化実験スクリプト
# 複数のパラメータ組み合わせで性能評価

# 環境設定
module load gcc/11.2.0
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_SCHEDULE=dynamic

echo "=== TSUBAME Dynamic Pruning Experiments ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK} cores"
echo "Start time: $(date)"
echo "============================================="

# 高性能版ビルド
echo "Building optimized version..."
make release
echo ""

# 実験パラメータ設定
ICB_DIR="ICB/ex01-new-4gen/"
SEEDS=(0 1 2)
BLOCK_LENGTHS=(3 5 10)
SPARSE_THRESHOLDS=("1e-5" "1e-6" "1e-7")

# 実験実行
experiment_count=0
total_experiments=$((${#SEEDS[@]} * ${#BLOCK_LENGTHS[@]} * ${#SPARSE_THRESHOLDS[@]}))

echo "Running ${total_experiments} experiments..."
echo ""

for seed in "${SEEDS[@]}"; do
    for block_length in "${BLOCK_LENGTHS[@]}"; do
        for threshold in "${SPARSE_THRESHOLDS[@]}"; do
            experiment_count=$((experiment_count + 1))

            echo "----------------------------------------"
            echo "Experiment ${experiment_count}/${total_experiments}"
            echo "Parameters: N=${block_length}, seed=${seed}, threshold=${threshold}"
            echo "Command: ./sim ${ICB_DIR} ${block_length} ${seed}"
            echo "Started at: $(date)"

            # 実行時間測定
            start_time=$(date +%s)
            ./sim ${ICB_DIR} ${block_length} ${seed} 2>&1 | while IFS= read -r line; do
                echo "[N=${block_length},seed=${seed},th=${threshold}] $line"
            done
            end_time=$(date +%s)

            execution_time=$((end_time - start_time))
            echo "Execution time: ${execution_time} seconds"
            echo "Completed at: $(date)"
            echo ""
        done
    done
done

echo "============================================="
echo "All experiments completed at: $(date)"
echo ""

echo "Summary:"
echo "  Total experiments: ${total_experiments}"
echo "  Block lengths tested: ${BLOCK_LENGTHS[*]}"
echo "  Seeds tested: ${SEEDS[*]}"
echo "  Thresholds tested: ${SPARSE_THRESHOLDS[*]}"
echo ""

echo "Analysis recommendations:"
echo "  1. Compare execution times across different N values"
echo "  2. Evaluate BER performance vs sparse threshold"
echo "  3. Check memory usage consistency (should be <100MB)"
echo "  4. Verify OpenMP scaling efficiency"
echo ""

echo "Results saved in: experiments_result_${SLURM_JOB_ID}.log"
echo "============================================="