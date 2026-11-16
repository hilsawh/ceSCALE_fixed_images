STEP 3 of pre-processing: Cellpose segmentation

Runs Cellpose in 3D on the DAPI TIFFs from Step 2 and saves one mask TIFF per embryo. It runs as a SLURM job on the cluster.

- **Runtime**: 28 seconds/embryo. The script requests 35 minutes, which covers about 70 embryos. Increase \--time if you have more.  
- **Requested:** gpus-per-node=1, nodes=1, cpus-per-task=4

&nbsp;

batch\_segment.sh is the script you edit and submit from the cluster with: sbatch batch\_segment.sh

The script loads python/3.11, gcc/12.3 and opencv/4.10.0, then activates a Python environment at \$HOME/cellpose. This environment must already exist with Cellpose installed.&nbsp;

&nbsp;

Before you submit, edit these:

- \--dir is the path to your "sampleTiffs/" folder from Step 2\.  
- \--savedir is the path where the masks are written, for example a "masks/" folder next to "sampleTiffs/".  
- /scratch/bazzi appears in the mkdir and export HOME line. Change it to your own scratch folder.  
- \--account is the lab's allocation on the cluster (e.g., def-YOUR\_PI).  
- \--anisotropy is scaled to each dataset's z-step. It is the z-step divided by the xy pixel size (110 nm).  
  - Dataset\_3 anisotropy: 150 / 110 \~= 1.4 (150 nm z-step)  
  - Dataset\_2 anisotropy: 200 / 110 \~= 1.8 (200 nm z-step)  
  - If your z-step or pixel size is different, calculate your own value.

&nbsp;

Output

One mask file per input TIFF is written to the \--savedir folder. Each mask is a Z stack where every nucleus has its own integer label and background is 0\.