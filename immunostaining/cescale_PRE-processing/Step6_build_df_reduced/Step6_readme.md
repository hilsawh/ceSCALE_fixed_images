STEP 6 of pre-processing: Build df\_reduced for CeSCALE

&nbsp;

Reads the text files from Step 5 and builds df\_reduced, a list with one table per embryo. This is the input for the CeSCALE notebooks 04\_a (1 P cell) and 04\_b (2 P cells).

&nbsp;

There are two notebooks. Use the one that matches your dataset:

- 1 P cell (dataset\_2): process\_ncoords\_df\_1Pcell.ipynb  
- 2 P cells (dataset\_3): process\_ncoords\_df\_2Pcell.ipynb

&nbsp;

**Runtime and resources:**

- process\_ncoords\_df\_1Pcell.ipynb  
  - Runtime: 2.4 minutes/embryo  
  - Requested: \-N 1 \--ntasks 8 \--cpus-per-task=2 \--mem 120G  
- process\_ncoords\_df\_2Pcell.ipynb  
  - Runtime: 5.8 minutes/embryo  
  - Requested: \-N 1 \--ntasks 8 \--cpus-per-task=2 \--mem 200G

&nbsp;

Open the notebook in Jupyter from the folder that contains "Ncoords\_1\_txt/" and run the cells from top to bottom. It needs numpy and pandas.

&nbsp;

**What the notebook does:**

- Reads one measurement table per embryo from "Ncoords\_1\_txt/". The 1 P cell notebook reads the nstats files. The 2 P cell notebook reads the nstats\_merged files, where the two P cells are one row at the end.  
- Keeps only the columns CeSCALE needs and renames the centroid columns to X, Y and Z.  
- Adds the P cell from the pcellnucID file (1 P cell) or the pcellMergednucID file (2 P cells).  
- Saves the list of tables as df\_reduced.

&nbsp;

**Output:**

- The last cell saves df\_reduced\_YOUR\_DATASET.pkl in the current folder. Change the file name for your dataset before you run it.  
- df\_reduced\_YOUR\_DATASET.pkl is the input for CeSCALE notebook 04\_a or 04\_b.  
- Move the file to CeSCALE's "segmentation\_embryos\_dataset/analysis/analysis\_files/segmentation\_results/" folder yourself. The file name that CeSCALE loads is written in the first cells of 04\_a and 04\_b, so change it there to match.