% scoring_CeSCALE_results_E.m
%
% Scores CeSCALE cell-identity assignments against lineage-tracing ground
% truth for one lineage. A cell is CeSCALE-positive if its assigned name
% starts with key_start_word (e.g. "E"). For each image, reports
% TP/FP/FN/TN, sensitivity and specificity in original Cellpose label IDs.
%
% Inputs
%   cescaleCSV : CeSCALE output, one row per cell. May hold several
%                replicates; sample_id looks like "R2_merged405_0_10_FOV1".
%                The script does not read a cell-ID column; it matches
%                cells by row order. Cellpose labels must be contiguous (1..N).
%                p_cells_num = 1: not for 1-P-cell CeSCALE output; use
%                scoring_CeSCALE_results_ABa.m for that
%                p_cells_num = 2: the 2 P-cells were merged into one cell and
%                placed in the last row before CeSCALE, so rows are
%                [other Cellpose labels in order, merged P-cell]. The script
%                reads which labels are P-cells from pcellCSV and converts
%                each row back to its Cellpose label(s).
%   lineageCSV : col 1 = mask/image name, col 2 = comma-separated Cellpose
%                labels that are lineage-positive (ground truth).
%   pcellCSV   : col 1 = mask/image name, col 2 = comma-separated Cellpose
%                labels of the P-cells.
%   Matching: image names in lineageCSV/pcellCSV (minus any "_cp_masks..."
%   suffix) are matched to sample_id (minus everything up to "merged").
%
% Output
%   output_file (.xlsx), one row per image.
%
% Some functions in this script were originally written with assistance 
% from ChatGPT (OpenAI) and were then checked manually.
%
% Requires MATLAB R2020a or later. Tested with R2025b

%% ---------------- Settings ----------------

% Positive class: CeSCALE names starting with this
key_start_word = "E";

% Number of P-cells. Keep at 2 (see cescaleCSV above).
p_cells_num = 2;

% Replicate to analyze. Use "" if sample_id has no replicate prefix
replicate_tag = "R3";

% Paths
data_dir    = '/path/to/data'; % point this to your data folder
cescaleCSV  = fullfile(data_dir, 'df_aligned_XYZ_complete.csv');
lineageCSV  = fullfile(data_dir, 'Lineage Numbers.csv');
pcellCSV    = fullfile(data_dir, 'P-cell Numbers.csv');
output_file = 'CeSCALE_Lineage_Classification_Results_E_R3.xlsx';

% CeSCALE column names
id_col         = 'sample_id';
label_col      = 'Assigned_name';
procrustes_col = 'min_Procrustes_mod_distance';

%% ---------------- Load CeSCALE output ----------------

T = readtable(cescaleCSV, 'VariableNamingRule', 'preserve');

% Keep only this replicate (must happen before grouping by image name)
T = T(startsWith(string(T.(id_col)), replicate_tag + "_"), :);
assert(height(T) > 0, 'No sample_id starts with "%s_". Check the prefix format.', replicate_tag);

assert(ismember(procrustes_col, T.Properties.VariableNames), ...
    'Column "%s" not found. Check the CSV header.', procrustes_col);

% "R2_merged405_0_10_FOV1" -> "405_0_10_FOV1"
cescale_file_names = regexprep(string(T.(id_col)), '^.*merged', '');
cescale_labels     = string(T.(label_col));
cescale_procrustes = T.(procrustes_col);

% Group by sample:
%   labelMap(id)      = CeSCALE names, in CeSCALE row order
%   procrustesMap(id) = Procrustes distance
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

    trueLabels    = lineageMap(key);   % ground-truth positives (Cellpose IDs)
    cescaleLabels = labelMap(key);     % CeSCALE names, in CeSCALE row order
    nLabels       = numel(cescaleLabels);
    predRows      = find(startsWith(cescaleLabels, key_start_word));

    % Convert CeSCALE-positive rows -> original Cellpose label IDs
    if p_cells_num == 2
        % The 2 P-cells were merged into one cell and placed in the last
        % row before CeSCALE. Rows 1..n-1 are the other Cellpose labels in
        % their original order; row n stands for both P-cells, e.g.
        %   original labels:   1 2 3 4 5 6   (P-cells = 2, 5)
        %   CeSCALE row order: 1 3 4 6 [2+5]
        allLabels = 1:(nLabels + 1);
        nonPcell  = setdiff(allLabels, pcellLabels, 'stable');
        if numel(pcellLabels) ~= 2 || numel(nonPcell) ~= nLabels - 1
            warning('P-cell mapping mismatch for key %s. Skipping.', key);
            continue
        end
        if any(trueLabels > nLabels + 1)
            warning('Lineage label above %d for key %s. Check the P-cell merge.', nLabels + 1, key);
        end
        predLabels = nonPcell(predRows(predRows < nLabels));
        if any(predRows == nLabels)
            predLabels = [predLabels, pcellLabels];
        end
    else
        % Row k = Cellpose label k
        allLabels  = 1:nLabels;
        predLabels = allLabels(predRows);
    end

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