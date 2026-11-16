STEP 4 of pre-processing: Mask correction in Napari

Opens each embryo in Napari with its Cellpose masks so you can fix the segmentation by hand, mark the P cells and the lineage-positive cells, and save the corrected masks.

&nbsp;

Run the script (Step4\_Mask\_Correction\_GUI.py) with python. It opens a Napari window in 3D with a "Mask correction" panel on the right. It needs napari, tifffile, numpy, scipy, scikit-image, natsort and qtpy.

&nbsp;

**Runtime and resources**: Run it on a local computer with a display. No cluster is needed. The work is done by hand, so the time depends on the user and on how many masks need to be fixed.

&nbsp;

**Loading data**

Use the buttons at the top of the panel to pick your folders:

- Select image folder: "sampleTiffs/" from Step 2\.  
- Select mask folder: "masks/" from Step 3\.  
- Select PG folder: "tiffs\_P\_Granules/" from Step 2\.  
- Select lineage fluorescence folder: "tiffs\_fluor/" from Step 2\.

Files are matched across folders by the sample and FOV numbers in their names, so DAPI\_0\_6\_FOV2.tif is paired with DAPI\_0\_6\_FOV2\_cp\_masks.tif. Images without a mask are left out.&nbsp;

&nbsp;

**Moving around**

- Left and right arrow keys go to the previous and next embryo.  
- "Switch to 2D" changes between the 3D view and single slices.  
- Click a mask to select it, click again to unselect. You can select several.  
- Type a mask number in "Find mask label" and press Enter to select it and jump to it.  
- C clears the selection. S shows only the DAPI signal inside the selected masks.  
- The layer list on the left of the Napari window lets you show or hide each layer (image, mask, P-Granules, lineage). Turning layers on and off makes it much easier to judge a mask.

&nbsp;

**Fixing masks**

- D deletes the selected masks.  
- M merges the selected masks into one (select two or more).  
- W splits one selected mask. Set "Number of Masks" first. A preview appears; click "Confirm watershed / new nuclei" to keep it. If the split is wrong, click "Enable point selection", place one point in each nucleus, and press W again.  
- "Min Distance" is the smallest distance in pixels allowed between the automatic split points. If it is too high, the split can give fewer masks than you asked for. Change "Min Distance" or "Number of Masks" and press W again.  
- \* N adds missing nuclei. Click "Enable point selection", place one point in the centre of each missing nucleus, then press N. "Missing Mask Threshold" is the DAPI intensity that counts as nucleus, and "Max Radius" limits how far each new mask can grow. Confirm the preview to keep it.  
- Place seed points in 2D view so they land where you intend. "Clear points" removes them.  
- You can re-run W or N with different points or settings as many times as you like before confirming.  
- &nbsp;"Delete Small Masks" removes every mask smaller than "Minimum Mask Size" (500 voxels by default).  
- "Dialate" grows the selected masks by "Pixel Dialation" pixels without overwriting their neighbours.  
- "Smooth all masks" rounds every mask while keeping its volume.  
- Ctrl+Z undoes the last edit, up to 10 edits.  
- Deleting or merging renumbers the masks, so mask numbers can change while you work.

&nbsp;

**P cells and lineage cells**

- When a PG folder is selected, P-granule signal is shown in magenta and the masks with the most P-granule signal are marked as P cells.&nbsp;  
- When a lineage folder is selected, the signal is shown in green and the masks with the most bright lineage signal are marked as lineage cells.&nbsp;  
- Set "Number of P-cells" (2 by default) and "Number of lineage cells" (8 by default) to match your embryo, and press Enter to run the detection again.  
  - **Note**: the number of P-cells per embryo for dataset\_3 was 2 (Z2, Z3) and for dataset\_2 was 1\. The box starts at 2 every time you open the GUI, so for dataset\_2 change it to 1 before you start.  
- The automatic choice is only a suggestion. It always picks the top-ranked masks up to the number you set, even if the embryo has fewer positive cells. Check every embryo, and fix mistakes by selecting a mask and clicking "Add selected to..." or "Remove selected from...".

&nbsp;

**Saving**

Ctrl+S saves the current embryo. Nothing is saved automatically. If you move to another embryo with unsaved edits, a warning appears; doing the same action a second time discards your edits.

&nbsp;

**Output**

- Everything is written to "Adjusted\_Masks/", created next to your mask folder.  
- One corrected mask file per embryo, named after the original mask file with \_adjusted added, for example DAPI\_0\_6\_FOV2\_cp\_masks\_adjusted.tif. Each nucleus has its own mask number. Numbers are 1 to N with no gaps.  
- "P-cell Numbers.csv" lists the mask numbers of the P cells for each embryo.  
- "Lineage Numbers.csv" lists the mask numbers of the lineage cells for each embryo.  
- Saving the same embryo again overwrites its mask and updates its row in both files. The original masks from Step 3 are never changed.

&nbsp;

**How the automatic steps work**

- **Splitting a mask**: if you placed seed points inside the selected mask, those are used. Otherwise the script picks the points furthest from the mask's edge, as many as "Number of Masks" and at least "Min Distance" apart. A watershed on the distance from the mask's edge then divides the mask between the points.  
- **Adding a missed mask**: a new mask can only use voxels that are brighter than "Missing Mask Threshold" and within "Max Radius" of a seed point. New masks do not overwrite existing ones.  
- **Automatic detection of P cells**: P-granule voxels below the 95th percentile of the image are set to zero. The mean P-granule signal is calculated for every mask, and the top 5 masks are kept. Each of these is grown by 5 pixels and the mean is calculated again, because P granules sit just outside the nucleus. The top masks, up to "Number of P-cells", are marked as P cells.  
- **Automatic detection of lineage cells**: the 90th percentile of lineage signal inside all masks is used as a threshold. For each mask, the fraction of its voxels above that threshold is calculated. The masks with the highest fractions, up to "Number of lineage cells", are marked as lineage cells.

&nbsp;

**Settings at the top of the script**

- MAX\_UNDO sets the number of undo steps.&nbsp;  
- The PERCENTILE settings control how much P-granule and lineage signal is displayed and counted.&nbsp;  
- LINEAGE\_DISPLAY\_MULTIPLIER brightens the lineage display.

&nbsp;

**Next step**

The corrected masks and "P-cell Numbers.csv" are the input for Step 5\. "Lineage Numbers.csv" is not used in Step 5\. It is used in the post-processing scoring script that calculates CeSCALE sensitivity and specificity.