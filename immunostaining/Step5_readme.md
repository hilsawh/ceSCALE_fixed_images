**STEP 5 of pre-processing: Nuclear coordinates (Ncoords) for CeSCALE**

Measures every nucleus in the corrected masks, marks which nuclei are P cells using "P-cell Numbers.csv" from Step 4, and writes the text files for Step 6\.

&nbsp;

There are two versions, one for embryos with 2 P cells and one for embryos with 1 P cell. Use the pair that matches your dataset. Each version has two scripts that you run in order from MATLAB:

- 2 P cells (dataset\_3): Step5a\_CeSCALE\_Ncoords\_2Pcell.m, then Step5b\_convert2csv\_2Pcell.m  
- 1 P cell (dataset\_2): Step5a\_CeSCALE\_Ncoords\_1Pcell.m, then Step5b\_convert2csv\_1Pcell.m

&nbsp;

**Runtime and resources:**

- This was run in an Open OnDemand MATLAB session on Trillium with 16 cores and 120 GB of memory requested.  
- dataset\_2 (1 P cell scripts): Step 5a took about 5 seconds per embryo. Step 5b took about 17 seconds per embryo.  
- dataset\_3 (2 P cell scripts): Step 5a took about 18 seconds per embryo. Step 5b took about 46 seconds per embryo.

&nbsp;

MATLAB's current folder must be the one that contains "Adjusted\_Masks/". The scripts need the Image Processing Toolbox. Both versions write to the same folder names, so keep each dataset in its own folder.

&nbsp;

**Step 5a: measure nuclei and mark P cells**

- Converts every mask TIFF in "Adjusted\_Masks/" to a .mat file in "matMasks/".  
- Measures every nucleus in 3D with regionprops3 (centroid, volume). Centroids are x, y, z in voxels.  
- Reads "P-cell Numbers.csv" and matches each row to its embryo by name, for example DAPI\_0\_6\_FOV2.  
- Saves one file per embryo in "Ncoords\_1/", for example DAPI\_0\_6\_FOV2.mat, and saves the main variables again as separate files in "Ncoords\_1/parsed/".

&nbsp;

What is saved in the 2 P cell version:

- nucL is the mask, nstats is the table of measurements, age is the number of nuclei.  
- z2cellnucID and z3cellnucID are the mask numbers of the two P cells. The lower mask number is called Z2 and the higher one Z3. This is only a naming rule; the script does not know which cell is really Z2 or Z3.  
- z2center and z3center are their centroids, and pcellncenter is the point halfway between them.  
- nstats\_separate is the full table with the two P cells as separate rows. somcellnucID\_separate lists the mask numbers of all other nuclei.  
- nstats\_merged is the table with the two P cells replaced by one row at the end, placed at pcellncenter. pcellMergednucID is the row number of that row, and somcellnucID\_merged lists the row numbers of all other nuclei. In the merged row only the centroid is meaningful. Volume, EquivDiameter and ConvexVolume are set to NaN.

&nbsp;

What is saved in the 1 P cell version:

- nucL, nstats and age as above.  
- pcellnucID is the mask number of the P cell and pcellncenter is its centroid.  
- somcellnucID lists the mask numbers of all other nuclei.

&nbsp;

**Step 5b: write text files for CeSCALE**

Reads every .mat file in "Ncoords\_1/" and writes its tables and ID lists as comma-separated .txt files in "Ncoords\_1\_txt/". Each file is named by the variable followed by the embryo name, for example nstats\_mergedDAPI\_0\_6\_FOV2.txt.

&nbsp;

- 2 P cell version writes nstats\_separate, nstats\_merged, pcellMergednucID, somcellnucID\_merged, somcellnucID\_separate, pcellncenter, z2cellnucID and z3cellnucID.  
- 1 P cell version writes nstats, pcellnucID, somcellnucID and pcellncenter.

&nbsp;

**Output**

"Ncoords\_1\_txt/" is the input for Step 6\. "matMasks/" and "Ncoords\_1/" are intermediate files.