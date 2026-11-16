#!/bin/bash

#SBATCH --time=0:35:00
#SBATCH --account=def-YOUR_PI
#SBATCH --gpus-per-node=1
#SBATCH --nodes=1
#SBATCH  --cpus-per-task=4


module load python/3.11 gcc/12.3 opencv/4.10.0
source $HOME/cellpose/bin/activate

mkdir -p /scratch/USERNAME/.cellpose && export HOME=/scratch/USERNAME

python -m cellpose --dir /scratch/USERNAME/YOUR_DATASET/sampleTiffs --anisotropy YOUR_ANISOTROPY --do_3D  --diameter 35 --verbose --pretrained_model cpsam --save_tif --use_gpu --no_npy --flow3D_smooth 0 --cellprob_threshold 0 --savedir /scratch/USERNAME/YOUR_DATASET/masks

