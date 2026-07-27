%% Post Cellpose/Napari CeSCALE Pre-processing for Embryos with 1 P-Cell

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
        continue
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

    % Convert "12" into 12
    pLabelStrings = strtrim(split(pLabelText, ","));
    pCellLabels = str2double(pLabelStrings);
    pCellLabels = pCellLabels(~isnan(pCellLabels));
    pCellLabels = pCellLabels(:);

    if numel(pCellLabels) ~= 1
        warning('Expected exactly 1 P-cell label for %s, found %d. Skipping.', ...
                imageBaseName, numel(pCellLabels));
        continue
    end

    pcellnucID = pCellLabels(1);

    % Make sure the P-cell label actually exists in the mask
    pcellRow = find(nucLabels == pcellnucID);
    
    if isempty(pcellRow)
        warning('P-cell label %d not found in nucLabels for %s. Skipping.', ...
                pcellnucID, imageBaseName);
        continue
    elseif numel(pcellRow) > 1
        warning('Multiple matches for P-cell label %d in %s. Skipping.', ...
                pcellnucID, imageBaseName);
        continue
    end

    pcellncenter = ncenters(pcellRow, :);
    somcellnucID = setdiff(allIDs, pcellnucID);

    %% Save all variables together

    outCombined = fullfile(statsOutDir, imageBaseName + ".mat");
    
    save(outCombined, ...
    'nucL', 'nstats', 'age', 'ncenters', 'nvols', ...
    'pcellnucID', 'somcellnucID', 'pcellncenter', ...
    'nucLabels', 'allIDs', ...
    'imageBaseName', ...
    '-v7.3');
    
    %% Save important variables separately
    
    save(fullfile(parsedOutDir, imageBaseName + "_nucL.mat"), 'nucL', '-v7.3');
    save(fullfile(parsedOutDir, imageBaseName + "_nstats.mat"), 'nstats');
    save(fullfile(parsedOutDir, imageBaseName + "_age.mat"), 'age');
    save(fullfile(parsedOutDir, imageBaseName + "_ncenters.mat"), 'ncenters');
    save(fullfile(parsedOutDir, imageBaseName + "_nvols.mat"), 'nvols');
    save(fullfile(parsedOutDir, imageBaseName + "_nucLabels.mat"), 'nucLabels');
    save(fullfile(parsedOutDir, imageBaseName + "_allIDs.mat"), 'allIDs');
    save(fullfile(parsedOutDir, imageBaseName + "_pcellnucID.mat"), 'pcellnucID');
    save(fullfile(parsedOutDir, imageBaseName + "_somcellnucID.mat"), 'somcellnucID');
    save(fullfile(parsedOutDir, imageBaseName + "_pcellncenter.mat"), 'pcellncenter');

    fprintf('Saved nuclear coordinate data for %s\n', imageBaseName);

end

clear pcellCSV maskMatDir statsOutDir parsedOutDir T csvImageNames csvPLabelsText