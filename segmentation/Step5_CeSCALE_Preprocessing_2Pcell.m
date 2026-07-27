%% Post Cellpose/Napari CeSCALE Pre-processing for Embryos with 2 P-Cells

maskDir = 'Adjusted_Masks';
maskfiles = dir(fullfile(maskDir, '*.tif'));
outMaskDir = 'matMasks';
if ~exist(outMaskDir, 'dir')
    mkdir(outMaskDir);
end

for k = 1:length(maskfiles)
    baseFileName = maskfiles(k).name;
    names = split(baseFileName,".");
    samplename = names(1,:);
    samplename = string(samplename);
    fullLoc = fullfile(maskDir, baseFileName);


    %here read in as an array of images
    imgMask = tiffreadVolume(fullLoc);
    filename = fullfile(outMaskDir, samplename + ".mat");
    save(filename, 'imgMask');

end
clear maskfiles

%% Step 2: Load P-cell CSV and calculate per-embryo nuclear stats

pcellCSV = fullfile(maskDir, 'P-cell Numbers.csv');   % <-- change if needed

maskMatDir = 'matMasks';
statsOutDir = 'Ncoords_1';
parsedOutDir = fullfile(statsOutDir, 'parsed');

if ~exist(statsOutDir, 'dir')
    mkdir(statsOutDir);
end

if ~exist(parsedOutDir, 'dir')
    mkdir(parsedOutDir);
end

% Read the whole file as raw text and parse manually
fid = fopen(pcellCSV, 'r');
lines = {};
while ~feof(fid)
    lines{end+1} = fgetl(fid);
end
fclose(fid);

csvImageNames = strings(length(lines)-1, 1);
csvPLabelsText = strings(length(lines)-1, 1);

for i = 2:length(lines)  % skip header
    % Find the first comma that separates image name from p-cell pair
    % Image names have underscores not commas, so first comma = split point
    firstComma = strfind(lines{i}, ',');
    if isempty(firstComma)
        continue
    end
    csvImageNames(i-1) = strtrim(lines{i}(1:firstComma(1)-1));
    csvPLabelsText(i-1) = strtrim(lines{i}(firstComma(1)+1:end));
end

csvImageNames = regexprep(csvImageNames, "_cp_masks.*", "");
csvPLabelsText = strrep(csvPLabelsText, '"', '');
disp(csvImageNames)
disp(csvPLabelsText)


% Loop through all MAT masks
maskMatFiles = dir(fullfile(maskMatDir, '*.mat'));

for k = 1:length(maskMatFiles)

    baseFileName = maskMatFiles(k).name;
    [~, maskBaseName, ~] = fileparts(baseFileName);

    % Example:
    % maskBaseName = "405_0_1_FOV1_cp_masks_adjusted"
    % imageBaseName = "405_0_1_FOV1"
    imageBaseName = regexprep(string(maskBaseName), "_cp_masks.*", "");

    fprintf('\nProcessing embryo: %s\n', imageBaseName);

    % Load mask
    maskFile = fullfile(maskMatDir, baseFileName);
    L = load(maskFile, 'imgMask');
    nucL = L.imgMask;

    % Regionprops3 nuclear stats
    nstats = regionprops3(nucL, 'all');

    % Remove invalid rows, if any
    if ismember('Centroid', nstats.Properties.VariableNames)
        nstats = nstats(~any(isnan(nstats.Centroid), 2), :);
    end

    if ismember('Volume', nstats.Properties.VariableNames)
        nstats = nstats(nstats.Volume > 0, :);
    end

    age = height(nstats);
    ncenters = nstats.Centroid;
    nvols = nstats.Volume;

    % Real label IDs in the mask
    nucLabels = unique(nucL(:));
    nucLabels(nucLabels == 0) = [];

    % Important check
    if numel(nucLabels) ~= height(nstats)
        warning(['Number of unique labels in nucL does not equal height(nstats) for %s. ' ...
                 'Check whether labels are non-consecutive or regionprops3 removed rows.'], ...
                 imageBaseName);
    end

    allIDs = nucLabels;

    %% Find matching P-cell row in CSV

    rowIdx = find(strcmpi(strtrim(csvImageNames), imageBaseName));

    if isempty(rowIdx)
        warning('No P-cell CSV entry found for %s. Skipping this embryo.', imageBaseName);
        continue
    elseif numel(rowIdx) > 1
        warning('Multiple P-cell CSV entries found for %s. Skipping this embryo.', imageBaseName);
        continue
    end

    pLabelText = csvPLabelsText(rowIdx);

    % Convert "12,13" into [12; 13]
    pLabelStrings = strtrim(split(pLabelText, ","));
    pCellLabels = str2double(pLabelStrings);
    pCellLabels = pCellLabels(~isnan(pCellLabels));
    pCellLabels = pCellLabels(:);

    if numel(pCellLabels) ~= 2
        warning('Expected exactly 2 P-cell labels for %s, found %d. Skipping.', ...
                imageBaseName, numel(pCellLabels));
        continue
    end

    z2cellnucID = pCellLabels(1);
    z3cellnucID = pCellLabels(2);

    % Match real P-cell labels to rows of nstats
    [tf, closestTwoIdx] = ismember(pCellLabels, nucLabels);

    if any(~tf)
        warning('Some P-cell labels for %s were not found in nucLabels. Skipping.', imageBaseName);
        disp(pCellLabels(~tf));
        continue
    end

    z2center = ncenters(closestTwoIdx(1), :);
    z3center = ncenters(closestTwoIdx(2), :);

    % Merged P-cell centroid
    pcellncenter = mean([z2center; z3center], 1);

    %% Build separate and merged versions

    % Somatic IDs when Z2 and Z3 are kept separate
    somcellnucID_separate = setdiff(allIDs, [z2cellnucID; z3cellnucID]);

    % Original stats table
    nstats_separate = nstats;

    % Merged stats table
    nstats_merged = nstats_separate;
    nstats_merged(closestTwoIdx, :) = [];

    % Use first P-cell row as template
    newRow = nstats_separate(closestTwoIdx(1), :);

    % Replace centroid with average P-cell center
    newRow.Centroid = pcellncenter;

    % These values are artificial after merging
    if ismember('Volume', nstats_merged.Properties.VariableNames)
        newRow.Volume = NaN;
    end

    if ismember('EquivDiameter', nstats_merged.Properties.VariableNames)
        newRow.EquivDiameter = NaN;
    end

    if ismember('ConvexVolume', nstats_merged.Properties.VariableNames)
        newRow.ConvexVolume = NaN;
    end

    % Append artificial merged P-cell row
    nstats_merged = [nstats_merged; newRow];

    % In merged table, P-cell is final row
    pcellMergednucID = height(nstats_merged);

    % Somatic IDs in merged table space
    somcellnucID_merged = setdiff((1:height(nstats_merged))', pcellMergednucID);

    %% Optional: split masks into P-cell and somatic masks

    pnucL = zeros(size(nucL), 'like', nucL);
    pnucL(nucL == z2cellnucID) = 1;
    pnucL(nucL == z3cellnucID) = 1;

    somnucL = nucL;
    somnucL(nucL == z2cellnucID) = 0;
    somnucL(nucL == z3cellnucID) = 0;

    %% Save all variables together

    outCombined = fullfile(statsOutDir, imageBaseName + ".mat");

    save(outCombined, ...
        'nucL', 'nstats', 'age', 'ncenters', 'nvols', ...
        'z2cellnucID', 'z3cellnucID', ...
        'pcellMergednucID', ...
        'somcellnucID_separate', 'somcellnucID_merged', ...
        'z2center', 'z3center', 'pcellncenter', ...
        'nstats_separate', 'nstats_merged', ...
        'nucLabels', 'allIDs', ...
        'pnucL', 'somnucL', ...
        'imageBaseName', ...
        '-v7.3');

    %% Save important variables separately

    save(fullfile(parsedOutDir, imageBaseName + "_nucL.mat"), 'nucL', '-v7.3');
    save(fullfile(parsedOutDir, imageBaseName + "_nstats.mat"), 'nstats');
    save(fullfile(parsedOutDir, imageBaseName + "_z2cellnucID.mat"), 'z2cellnucID');
    save(fullfile(parsedOutDir, imageBaseName + "_z3cellnucID.mat"), 'z3cellnucID');
    save(fullfile(parsedOutDir, imageBaseName + "_pcellMergednucID.mat"), 'pcellMergednucID');
    save(fullfile(parsedOutDir, imageBaseName + "_somcellnucID_separate.mat"), 'somcellnucID_separate');
    save(fullfile(parsedOutDir, imageBaseName + "_somcellnucID_merged.mat"), 'somcellnucID_merged');
    save(fullfile(parsedOutDir, imageBaseName + "_z2center.mat"), 'z2center');
    save(fullfile(parsedOutDir, imageBaseName + "_z3center.mat"), 'z3center');
    save(fullfile(parsedOutDir, imageBaseName + "_pcellncenter.mat"), 'pcellncenter');
    save(fullfile(parsedOutDir, imageBaseName + "_nstats_separate.mat"), 'nstats_separate');
    save(fullfile(parsedOutDir, imageBaseName + "_nstats_merged.mat"), 'nstats_merged');

    fprintf('Saved nuclear coordinate data for %s\n', imageBaseName);

end

clear pcellCSV maskMatDir statsOutDir parsedOutDir T csvImageNames csvPLabelsText