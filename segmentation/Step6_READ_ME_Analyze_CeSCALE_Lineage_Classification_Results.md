# Purpose of Code
This MATLAB file produces a .xlsx file with the sensitivity and specificity of CeSCALE cell classification relative to the fluorescently labeled cell ground truth. This helps to evaluate how effective CeSCALE is at identifying cell types in fixed embryos. 

# How to Run Code
1. First the user must adjust the key_start_word variable so it reflects a unique starting string of the cell type that they have fluorescently labeled. For example the starting string "E" is unique to intestinal cells. 
2. The user must adjust the number of P-cells in their embryos. Must either be 1 or 2.
3. The user must input the full path to the P-cell CSV, the lineage tracing CSV, and the CSV produced as output to CeSCALE.

Note: these steps should be enough to run the code correctly but sometimes the naming does not match between the CeSCALE CSV and the P-cell and lineage cell CSVs. In this case additional alterations may be needed to have the keys match. 
