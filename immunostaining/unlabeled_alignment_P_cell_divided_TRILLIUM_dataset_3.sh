#!/bin/bash
#SBATCH --job-name=unlabeled_alignment_P_cell_divided
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --output=logs/unlabeled_data_alignment_P_cell_divided.%j.o
#SBATCH --error=logs/unlabeled_data_alignment_P_cell_divided.%j.e
#SBATCH --mail-type=END,FAIL,TIME_LIMIT

if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
	cd "${SLURM_SUBMIT_DIR}"
else
	HERE="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
	cd "$(dirname "$HERE")"
fi
[[ -f scripts/unlabeled_data_alignment_P_cell_divided.py ]] || { echo "ERROR: not in the repo root ($PWD)" >&2; exit 1; }
mkdir -p logs

module load python/3.10
source /home/bazzi/cescale_latest_env/bin/activate || exit 1

export PYTHONPATH="$PWD${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=4
export MPLCONFIGDIR=/scratch/bazzi/.config/matplotlib

if [[ -f .env ]]; then
	source .env
fi
data_dir="${Cell_lineage_segmentation_DATA}"
save_dir="${Cell_lineage_PATH}"

if [[ -z "${data_dir}" || -z "${save_dir}" ]]; then
	echo "ERROR: set Cell_lineage_segmentation_DATA and Cell_lineage_PATH in $PWD/.env" >&2
	exit 1
fi

python -c "import sys, numpy, pandas, scipy, ot, cl_id; print('python', sys.executable, '| numpy', numpy.__version__, '| pandas', pandas.__version__, '| scipy', scipy.__version__, '| POT', ot.__version__, '| cl_id', list(cl_id.__path__))" || exit 1

python -u scripts/unlabeled_data_alignment_P_cell_divided.py \
	--data_dir "${data_dir}" --save_dir "${save_dir}" \
	--test-file df_reduced.pkl