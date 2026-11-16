**POST-PROCESSING: CeSCALE scoring and confidence overlays**

&nbsp;

These scripts are run after CeSCALE.&nbsp;

The scoring scripts compare CeSCALE's cell names with the lineage-positive cells you marked in Step 4 of pre-processing, and report sensitivity and specificity for each embryo.&nbsp;

The overlay scripts draw CeSCALE's cell names and per-cell confidence scores on the images of one embryo.

&nbsp;

**Runtime and resources:**

- The scripts do not need a cluster and can be run on a local computer.  
- Each script takes less than 1 minute per run.&nbsp;  
- A scoring run covers all embryos of one replicate.&nbsp;  
- An overlay run makes the figure panels for one embryo.  
- build\_viridis\_csv.m is run once.

&nbsp;

Use the scripts that match your dataset:

- 1 P cell, ABa lineage marker (dataset\_2): **scoring\_CeSCALE\_results\_ABa.m** and **confidence\_overlay\_ABa.m**  
- 2 P cells, E lineage marker (dataset\_3): **scoring\_CeSCALE\_results\_E.m** and **confidence\_overlay\_E.m**  
- For both: **build\_viridis\_csv.m**, run once before your first overlay

&nbsp;

All five are MATLAB scripts. Open the script, edit the Settings section at the top, and run it. Scoring and overlays do not depend on each other, so you can run them in either order.

&nbsp;

Row order and mask numbers:

- 1 P cell: row k is mask number k in "Adjusted\_Masks/".  
- 2 P cells: row k is not mask number k. The two P cells were merged into one row at the end in Step 5a, so every nucleus after the first P cell moves up. The E script reads the two P cell mask numbers from "P-cell Numbers.csv" and converts each row back to its mask number. For example, with masks 1 2 3 4 5 6 and P cells 2 and 5, CeSCALE's rows are 1 3 4 6 and then the merged P cell.  
- Masks must be numbered 1 to N with no gaps.  
- Do not sort or filter the rows of an embryo in the CeSCALE file. A changed row order gives wrong scores.

&nbsp;

**SCORING: scoring\_CeSCALE\_results\_ABa.m and scoring\_CeSCALE\_results\_E.m**

&nbsp;

What the scripts do for each embryo:

- A cell is CeSCALE-positive if its Assigned\_name starts with key\_start\_word, "ABa" or "E". This includes all daughters, for example ABal and ABar, or Ea and Ep.  
- A cell is lineage-positive if its mask number is listed for that embryo in "Lineage Numbers.csv".  
- TP is positive in both. FP is CeSCALE-positive only. FN is lineage-positive only. TN is positive in neither.  
- Sensitivity is TP / (TP \+ FN). Specificity is TN / (TN \+ FP).  
- In the E script, Z2 and Z3 count as two cells, although CeSCALE saw them as one.

&nbsp;

Inputs:

- df\_aligned\_XYZ\_complete.csv is the final output of CeSCALE notebook 04\_a (1 P cell) or 04\_b (2 P cells), saved in "analysis\_files/segmentation\_results/alignment\_results/". It has one row per cell. The scripts use the columns sample\_id, Assigned\_name and min\_Procrustes\_mod\_distance.  
- "Lineage Numbers.csv" from Step 4\. Column 1 is the mask file name. Column 2 is the mask numbers of the lineage-positive cells.  
- "P-cell Numbers.csv" from Step 4\. Column 1 is the mask file name. Column 2 is the mask number(s) of the P cell(s). The ABa script only copies these into the output. The E script needs them to convert rows back to mask numbers.

&nbsp;

Matching embryo names

Mask file names lose everything from "\_cp\_masks" onwards, so DAPI\_0\_6\_FOV2\_cp\_masks.tif becomes DAPI\_0\_6\_FOV2. This must be the same as sample\_id in the CeSCALE file after the script has cleaned it. Embryo names can repeat between dataset replicates. In the E script, run one replicate at a time, each with its own "Lineage Numbers.csv", "P-cell Numbers.csv" and output file name.

&nbsp;

Settings to edit

- key\_start\_word is the start of the cell names that count as positive.  
- data\_dir, cescaleCSV, lineageCSV and pcellCSV are the folder and the three input files.  
- output\_file is the name of the generated Excel file. It is written to MATLAB's current folder.

&nbsp;

**Output**

One Excel file with one row per embryo:

- ImageName is the embryo name.  
- PcellLabels is the mask number(s) of the P cell(s).  
- LineagePositiveLabels is the mask numbers of the lineage-positive cells.  
- CeSCALEPositiveLabels\_AdjustedToOriginalIDs is the mask numbers of the CeSCALE-positive cells.  
- Procrustes is min\_Procrustes\_mod\_distance. It is one value per embryo.  
- TP, FP, FN, TN, Sensitivity and Specificity.

&nbsp;

&nbsp;

**CONFIDENCE OVERLAYS: build\_viridis\_csv.m, confidence\_overlay\_ABa.m and confidence\_overlay\_E.m**

&nbsp;

build\_viridis\_csv.m

Run this once. It downloads matplotlib's colour map file (\_cm\_listed.py, version v3.11.2) from GitHub, reads the 256 viridis colours and saves them as viridis256.csv in MATLAB's current folder. The file has 256 rows of R, G, B values between 0 and 1\.

&nbsp;

THE OVERLAY SCRIPTS

Each run makes the figure panels for one embryo. viridis256.csv must be in MATLAB's current folder. Both overlay scripts need MATLAB R2020a or later. The ABa script also needs the Image Processing Toolbox.

&nbsp;

Inputs for **confidence\_overlay\_ABa.m**

- Everything is under dataRoot:  
- combined\_per\_cell\_confidence.csv is the per-cell confidence table from CeSCALE notebooks 06\_a or 06\_b (1Pcell or 2Pcell).  
- "Ncoords\_1/parsed/SampleNFOVncoords.mat" is the crop box from Step 2 in pre-processing.  
- "Cropped\_Embryo\_Images/DAPI\_0\_N\_FOVn.mat" is the cropped DAPI stack from Step 2 in pre-processing.&nbsp;  
- "hyb0\_channels/PGranule\_0\_N.mat" and "hyb0\_channels/Lineage\_0\_N.mat" are the full stacks from Step 1 in pre-processing.

&nbsp;

Inputs for **confidence\_overlay\_E.m**

- The confidence table is under dataRoot and its name is set in confFileName.&nbsp;  
- The images are under imageRoot:  
  - "Ncoords\_1/DAPI\_0\_N\_FOVn.mat" is the Step 5a pre-processing file for the embryo. The script reads z2center and z3center from it.  
  - "Ncoords\_1/parsed/SampleNFOVncoords.mat" is the crop box from Step 2 in pre-processing.  
  - "hyb0\_channels/DAPI\_0\_N.mat", "hyb0\_channels/PGranule\_0\_N.mat" and "hyb0\_channels/Lineage\_0\_N.mat" are the full stacks from Step 1 in pre-processing.

&nbsp;

&nbsp;

Settings to edit:

- dataRoot (and imageRoot and confFileName in the E script), runTag, SampleNum and FOV.  
- pixelSize is in micrometres per pixel (0.11) and scaleBarLength\_um is the length of the scale bar (5).  
- confRange is the range of the colour bar, \[0 0.7\]. Scores are coloured on a linear scale with viridis. Scores of 0.7 or more all get the top colour, and the last tick of the colour bar reads "≥0.7".  
- dotSize and the colours of the markers.  
- ABa script only: medK is the size of the median filter, \[3 3 3\].  
- E script only: eSlice and eHalfWidth. The lineage panel is a maximum projection of the planes from eSlice \- eHalfWidth to eSlice \+ eHalfWidth (set as 125 ± 10 in the script). Choose eSlice for each embryo.

&nbsp;

How the images are made

- ABa script: each stack gets a 3 by 3 by 3 median filter, then a maximum projection over all planes. The DAPI panel uses the stack from Step 2 of pre-processing. The P granule and lineage panels use the stacks from Step 1 of pre-processing, cropped to the crop box.  
- E script: no filter is applied. The stacks of this dataset were denoised with Denoise.ai in NIS-Elements before Step 1 of pre-processing. DAPI and P granule are maximum projections over all planes, cropped to the crop box. The lineage panel uses only the planes around eSlice.

&nbsp;

**Panels**

- DAPI with confidence. Every nucleus is a dot at its centroid (Initial\_X, Initial\_Y), coloured by its confidence score. Germline cells are not scored. The scale bar is at the bottom right. The title shows the embryo name and its Procrustes value.  
- P granules. The ABa script marks the P cell with a green diamond. The E script marks Z2 and Z3 with green diamonds labelled "Z2/3".  
- Lineage marker. Cells whose Assigned\_name starts with "ABa" (or "E") are marked with dots and their names. The ABa title gives the number of ABa cells. The E title gives the planes that were projected.  
- Colour bar, saved as its own file.

&nbsp;

**Output**

Files are written to a folder (overlays/) in MATLAB's current folder, which is made if it does not exist. Each panel is saved as a .tif and as a vector .pdf with a transparent background.