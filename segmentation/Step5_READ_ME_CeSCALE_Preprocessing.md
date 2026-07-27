# Purpose of Code
Calculates Nstats for each embryo's set of nuclei masks. Importantly this includes the centroid of each nucleus, which will be used for CeSCALE. There are 2 different versions of the code: one for embryos with 1 P-cell and one for embryos with 2 P-cells. If your embryos have a combination of 1 and 2 P-cells it is important that they are split into 2 separate folders, this can easily be done using the Napari GUI, and run separately. 

# How to Run Code
On line 3 change the maskDir to point to your "Adjusted_Masks" directory made in Napari. This folder should contain both the adjusted masks as well as the CSV file for P-cell numbers (exact name: 'P-cell Numbers.csv'). 
