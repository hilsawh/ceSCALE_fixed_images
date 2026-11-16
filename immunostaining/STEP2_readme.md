STEP 2 of pre-processing: Embryo outlining and cropping

Lets you outline each embryo by hand on the DAPI image, then crops every channel to each embryo and saves TIFF stacks for the downstream Cellpose and Napari steps (STEPS 3 and 4, respectively).

&nbsp;

**Runtime and resources:**

- This was run in an Open OnDemand MATLAB session on Trillium with 16 cores and 120 GB of memory requested.&nbsp;  
- Block 1 of the script is done by hand, so its time depends on the user.  
- dataset\_2: Blocks 2 to 4 took about 17 seconds per embryo.  
- dataset\_3: Blocks 2 to 4 took about 49 seconds per embryo.

&nbsp;

**Step2\_Segment\_Embryos.mlx** is a MATLAB live script with four blocks. Run them in order. It reads the files made in Step 1, so MATLAB's current folder must be the one that contains "hyb0\_channels/".&nbsp;

- NOTE: this script needs the Image Processing Toolbox.

&nbsp;

*Block 1: Outline embryos*

Opens the DAPI maximum projection for each sample in a pop-up window, one sample at a time. Click "Draw embryo", draw an ellipse around one embryo, adjust it, then double-click it. Click "Accept ellipse" to save it or "Redraw" to try again. Click "Draw another" for the next embryo in the same image, or "Next sample" to move on.

For each embryo, SampleNFOVncoords.mat holds the crop box (variables x1, x2, y1, y2).

&nbsp;

*Block 2: DAPI for Cellpose and Napari*

Crops the DAPI stack to each individual embryo chosen in Block 1, sets everything outside the ellipse to zero, subtracts background (disk radius 30), and smooths with a 5 by 5 by 5 median filter. Saves the result as DAPI\_0\_Sample\_FOVn.mat in "Cropped\_Embryo\_Images/" and as DAPI\_0\_Sample\_FOVn.tif in "sampleTiffs/". The TIFFs are the input for Cellpose (STEP 3), and are also one of the inputs for Napari (STEP 4).

&nbsp;

*Block 3: P-granule channel for Napari*

Crops the P-granule stack to each embryo's crop box. Saves PGranule\_0\_Sample\_FOVn.tif in "tiffs\_P\_Granules/". The TIFFs are one of the inputs for Napari (STEP 4).

&nbsp;

*Block 4: Lineage channel for Napari*

Crops the lineage stack to each embryo's crop box and subtracts background (disk radius 30). Saves Lineage\_0\_Sample\_FOVn.tif in "tiffs\_fluor/". The TIFFs are one of the inputs for Napari (STEP 4).

&nbsp;