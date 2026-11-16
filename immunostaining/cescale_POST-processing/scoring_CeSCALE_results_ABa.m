% scoring_CeSCALE_results_ABa.m
%
% Scores CeSCALE cell-identity assignments against lineage-tracing ground
% truth for one lineage. A cell is CeSCALE-positive if its assigned name
% starts with key_start_word (e.g. "ABa"). For each image, reports
% TP/FP/FN/TN, sensitivity and specificity in original Cellpose label IDs.
%
% Inputs
%   cescaleCSV : CeSCALE output (df_aligned_XYZ_complete.csv), one row per
%                cell. The script does not read a cell-ID column; it matches
%                cells by row order: row k of each sample_id = Cellpose
%                label k. Cellpose labels must be contiguous (1..N).
%                Assumes one P-cell per embryo, left in its original
%                position.
%   lineageCSV : col 1 = mask/image name, col 2 = comma-separated Cellpose
%                labels that are lineage-positive (ground truth).
%   pcellCSV   : col 1 = mask/image name, col 2 = Cellpose label of the
%                P-cell. Only reported in the output (PcellLabels column);
%                not used for scoring.
%   Image names in lineageCSV/pcellCSV are matched to sample_id after
%   removing any "_cp_masks..." suffix.
%
% Output
%   output_file (.xlsx), one row per image.
%
% Some functions in this script were originally written with assistance 
% from ChatGPT (OpenAI) and were then checked manually.
%
% Requires MATLAB R2017a or later. Tested with MATLAB R2025b

%% ---------------- Settings ----------------

% Positive class: CeSCALE names starting with this
key_start_word = "ABa";

% Paths
data_dir    = '/path/to/data'; % point this to your data folder
cescaleCSV  = fullfile(data_dir, 'df_aligned_XYZ_complete.csv');
lineageCSV  = fullfile(data_dir, 'Lineage Numbers.csv');
pcellCSV    = fullfile(data_dir, 'P-cell Numbers.csv');
output_file = 'CeSCALE_Lineage_Classification_Results.xlsx';

% Optional run/contributor prefix on sample_id.
% "" keeps all rows unchanged.
cescale_prefix = "";

% CeSCALE column names
id_col         = 'sample_id';
label_col      = 'Assigned_name';
procrustes_col = 'min_Procrustes_mod_distance';

%% ---------------- Load CeSCALE output ----------------

T = readtable(cescaleCSV);
T = T(startsWith(string(T.(id_col)), cescale_prefix), :);

cescale_file_names = erase(string(T.(id_col)), cescale_prefix);
cescale_labels     = string(T.(label_col));
cescale_procrustes = T.(procrustes_col);

% Group by sample:
%   labelMap(id)      = CeSCALE names, in row order
%   procrustesMap(id) = mod Procrustes distance (one value per sample; taken
%                       from its first row)
labelMap      = containers.Map('KeyType', 'char', 'ValueType', 'any');
procrustesMap = containers.Map('KeyType', 'char', 'ValueType', 'any');

for i = 1:numel(cescale_file_names)
    key = char(cescale_file_names(i));
    if isKey(labelMap, key)
        labelMap(key) = [labelMap(key), cescale_labels(i)];
    else
        labelMap(key)      = cescale_labels(i);
        procrustesMap(key) = cescale_procrustes(i);
    end
end

%% ---------------- Load ground truth ----------------

lineageMap = makeLabelMapFromCSV(lineageCSV);
pcellMap   = makeLabelMapFromCSV(pcellCSV);

%% ---------------- Score each image ----------------

outRows     = {};
lineageKeys = keys(lineageMap);

for k = 1:numel(lineageKeys)

    key = lineageKeys{k};

    if ~isKey(labelMap, key)
        warning('Key %s found in lineageMap but not in labelMap. Skipping.', key);
        continue
    end

    if isKey(pcellMap, key)
        pcellLabels = pcellMap(key);
    else
        warning('Key %s found in lineageMap but not in pcellMap.', key);
        pcellLabels = [];
    end

    trueLabels    = lineageMap(key);   % ground-truth positives
    cescaleLabels = labelMap(key);     % CeSCALE names; row k = Cellpose label k
    allLabels     = 1:numel(cescaleLabels);

    % CeSCALE-positive cells, as Cellpose label IDs
    predRows   = find(startsWith(cescaleLabels, key_start_word));
    predLabels = allLabels(predRows);

    % Confusion counts
    TP = numel(intersect(predLabels, trueLabels));
    FP = numel(setdiff(predLabels, trueLabels));
    FN = numel(setdiff(trueLabels, predLabels));
    TN = numel(setdiff(allLabels, union(predLabels, trueLabels)));

    % Metrics (0/0 returns NaN when a class is empty)
    sensitivity = TP / (TP + FN);
    specificity = TN / (TN + FP);

    outRows(end+1, :) = { ...
        string(key), ...
        listToText(pcellLabels), ...
        listToText(trueLabels), ...
        listToText(predLabels), ...
        procrustesMap(key), ...
        TP, FP, FN, TN, ...
        sensitivity, specificity};
end

if isempty(outRows)
    error('No images matched. Lineage keys: %s | CeSCALE keys: %s', ...
        strjoin(lineageKeys, ', '), strjoin(keys(labelMap), ', '));
end

resultsTable = cell2table(outRows, 'VariableNames', { ...
    'ImageName', 'PcellLabels', 'LineagePositiveLabels', ...
    'CeSCALEPositiveLabels_AdjustedToOriginalIDs', 'Procrustes', ...
    'TP', 'FP', 'FN', 'TN', 'Sensitivity', 'Specificity'});

writetable(resultsTable, output_file, 'FileType', 'spreadsheet');

%% ---------------- Helper functions ----------------

function txt = listToText(v)
    % Numeric vector -> "1,2,3" ("" if empty)
    if isempty(v)
        txt = "";
    else
        txt = join(string(v), ',');
    end
end

function labelMap = makeLabelMapFromCSV(inputFile)
    % Reads a CSV where:
    %   column 1 = image/mask filename
    %   column 2 = comma-separated label list, e.g. "12,13"
    % Returns a map: labelMap('405_0_11_FOV1') = [12 13]

    fid = fopen(inputFile, 'r');
    if fid == -1
        error('Could not open file: %s', inputFile);
    end

    lines = {};
    while ~feof(fid)
        lines{end+1} = fgetl(fid);
    end
    fclose(fid);

    csvImageNames = strings(length(lines)-1, 1);
    csvLabelText  = strings(length(lines)-1, 1);

    % Parse lines, skipping header; split on the first comma only
    for i = 2:length(lines)
        firstComma = strfind(lines{i}, ',');
        if isempty(firstComma)
            continue
        end
        csvImageNames(i-1) = strtrim(lines{i}(1:firstComma(1)-1));
        csvLabelText(i-1)  = strtrim(lines{i}(firstComma(1)+1:end));
    end

    % Strip Cellpose mask suffix and CSV quotes
    csvImageNames = regexprep(csvImageNames, "_cp_masks.*", "");
    csvLabelText  = strrep(csvLabelText, '"', '');

    labelMap = containers.Map('KeyType', 'char', 'ValueType', 'any');

    for i = 1:numel(csvImageNames)
        if strlength(csvImageNames(i)) == 0 || strlength(csvLabelText(i)) == 0
            continue
        end
        key = char(csvImageNames(i));
        labels = str2double(split(csvLabelText(i), ','));
        labels = labels(~isnan(labels));
        labelMap(key) = labels(:)';
    end
end