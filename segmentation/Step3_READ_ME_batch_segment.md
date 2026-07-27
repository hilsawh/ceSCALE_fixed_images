# Purpose of Code
This code performs batch nuclei segmentation of all DAPI images in a folder using Cellpose-SAM. These masks are then used downstream to determine the centroids of nuclei for CeSCALE. This code has been designed to run on a computational cluster, specifically Alliance Canada Trillium.

# How to Run Code
1. Change the path on line 11 to reflect the path to your own cellpose/bin/activate folder
2. Change the Home= variable on line 13 to correspond to the path that contains your cellpose folder
3. Change the MPLCONFIGDIR= variable on line 14 to correspond to the path to your own .config/matplotlib folder
4. Change the path immediately following --dir on line 16 to the folder containing your DAPI images saved as TIFFs in the previous step, typically saved as a folder called "sampleTiffs".
5. Change the path immediately following --savedir on line 16 to an existing folder where you would like to save your cellpose masks. 
