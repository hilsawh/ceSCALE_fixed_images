# Purpose of Code
### Code Block 1: Embryo Outline Using DAPI Channel
In code block 1 the user views the DAPI channel and draws ovals around the embryos they would like to use for analysis. These ovals are then used to mask background signal and crop the image around the embryos of interest.
### Code Block 2: DAPI preprocessing and Cellpose TIFF creation
In code block 2 the DAPI channel is cropped and strel background subtraction and median (5,5,5) smoothing is performed. The cropped and processed DAPI channel images are saved as .mat and .tiff files (used for Cellpose).
### Code Block 3: P-Granule Channel Processing and TIFF Creation for Napari
In code block 3 the P-Granule channel is cropped and converted to .tiff files so that they can be used in the Napari GUI to determine which nuclei masks correspond to the P-cells.
### Code Block 4: Lineage Channel Processing and TIFF Creation for Napari
In code block 4 the lineage fluorescence channel is cropped, strel background subtracted, and converted to .tiff files. These will also be used in the Napari GUI to determine which nuclei masks correspond to lineage fluorescent cells. 

# How to Run Code
Run each code section in order using the "Run Section" button in MATLAB. You should not have to change any variables within the code in any of the blocks. In the first block there will be a pop up where you will be able to view the max projection of the DAPI signal. You then click "Draw Embryo" and outline the embryo using an oval. Once you do the initial outline you can adjust it to better fit the embryo. When you are happy with the outline you can double click the embryo. If there is another embryo in the sample that you would like to outline, you click "Draw Another". If you are done with the sample click "Next Sample". Once you have completed the outlining for all embryos the first block of code will end. The 2nd, 3rd, and 4th blocks of code you just need to click "Run Section". The 3rd and 4th blocks of code are optional and should only be run if you have P-Granule and Lineage Fluorescence channels respectively. 
