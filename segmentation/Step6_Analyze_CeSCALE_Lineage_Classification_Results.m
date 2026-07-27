% Cell Type Key Start Word
key_start_word = "E";

% Number of P-Cells, either 1 or 2, if 2 P-cell CSV used to reorginize the
% CeSCALE IDs to match the Napari IDs
p_cells_num = 2;

% Lineage CSV
lineageCSV = 'C:\Users\james\OneDrive\Documents\Lineage Classification Test\Adjusted_Masks\Lineage Numbers.csv';
% CeScale CSV
T = readtable('C:\Users\james\OneDrive\Documents\Lineage Classification Test\df_aligned_XYZ_complete.csv');
% P-Cell CSV
pcellCSV = 'C:\Users\james\OneDrive\Documents\Lineage Classification Test\Adjusted_Masks\P-cell Numbers.csv';

% Creates dictionaries with keys being the image name.
% labelMap stores the list of CeSCALE labels ordered by Cellpose label ID.
% procrustesMap stores one Procrustes value per embryo.
cescale_file_names = string(T{:, 2});
cescale_file_names = erase(cescale_file_names, "_merged");

cescale_labels = string(T{:, 15});

% Procrustes distance measure from column 25
cescale_procrustes = T{:, 25};

labelMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
procrustesMap = containers.Map('KeyType', 'char', 'ValueType', 'any');

for i = 1:numel(cescale_file_names)

    key = char(cescale_file_names(i));
    label = cescale_labels(i);

    if isKey(labelMap, key)
        labelMap(key) = [labelMap(key), label];
    else
        labelMap(key) = label;
        procrustesMap(key) = cescale_procrustes(i);
    end

end

lineageMap = makeLabelMapFromCSV(lineageCSV);
pcellMap = makeLabelMapFromCSV(pcellCSV);

% Output rows
outRows = {};

lineageKeys = keys(lineageMap);

for k = 1:numel(lineageKeys)

    key = lineageKeys{k};

    % Only analyze files present in both dictionaries
    if ~isKey(labelMap, key)
        warning('Key %s found in lineageMap but not in labelMap. Skipping.', key);
        continue
    end

    % Get P-cell labels for this image
    if isKey(pcellMap, key)
        pcellLabels = pcellMap(key);
        pcellLabels = pcellLabels(:)';   % force row vector
    else
        warning('Key %s found in lineageMap but not in pcellMap. Assuming no P-cells.', key);
        pcellLabels = [];
    end

    % Ground truth positive labels from lineage tracing
    trueLabels = lineageMap(key);
    trueLabels = trueLabels(:)';

    % CeSCALE cell-type labels ordered by CeSCALE-adjusted label number
    cescaleLabels = string(labelMap(key));

    % ------------------------------------------------------------
    % Build mapping from CeSCALE index back to original Cellpose label ID
    %
    % Assumption:
    % P-cells were removed from their original positions and appended
    % to the end before CeSCALE.
    %
    % Example:
    % original labels:        1 2 3 4 5
    % pcellLabels:           3
    % CeSCALE index order:   1 2 4 5 3
    %
    % Therefore:
    % cescaleIndexToOriginalLabel(3) = 4
    % ------------------------------------------------------------

    if p_cells_num == 2
        nLabels = numel(cescaleLabels) + 1;
        allOriginalLabels = 1:nLabels;
        nonPcellOriginalLabels = setdiff(allOriginalLabels, pcellLabels, 'stable');
        cescaleIndexToOriginalLabel = [nonPcellOriginalLabels, pcellLabels];
    else
        nLabels = numel(cescaleLabels);
        cescaleIndexToOriginalLabel = 1:nLabels;
    end

    % Safety check
    if numel(cescaleIndexToOriginalLabel) ~= nLabels
        warning('Mapping length mismatch for key %s. Skipping.', key);
        continue
    end

    % CeSCALE predicted positive indices in CeSCALE-adjusted ordering
    predCescaleIndices = find(startsWith(cescaleLabels, key_start_word));

    % Convert CeSCALE indices back to original Cellpose/lineage label IDs
    predLabels = cescaleIndexToOriginalLabel(predCescaleIndices);

    % All labels considered, in original Cellpose/lineage label ID space
    allLabels = cescaleIndexToOriginalLabel;

    % Confusion counts using original label IDs
    TP = numel(intersect(predLabels, trueLabels));
    FP = numel(setdiff(predLabels, trueLabels));
    FN = numel(setdiff(trueLabels, predLabels));
    TN = numel(setdiff(allLabels, union(predLabels, trueLabels)));

    % Metrics
    if TP + FN == 0
        sensitivity = NaN;
    else
        sensitivity = TP / (TP + FN);
    end

    if TN + FP == 0
        specificity = NaN;
    else
        specificity = TN / (TN + FP);
    end

    % Convert positive-label lists to CSV-friendly strings
    if isempty(trueLabels)
        trueLabelsText = "";
    else
        trueLabelsText = join(string(trueLabels), ',');
    end

    if isempty(predLabels)
        predLabelsText = "";
    else
        predLabelsText = join(string(predLabels), ',');
    end

    if isempty(pcellLabels)
        pcellLabelsText = "";
    else
        pcellLabelsText = join(string(pcellLabels), ',');
    end

    % Add row
    outRows(end+1, :) = { ...
        string(key), ...
        pcellLabelsText, ...
        trueLabelsText, ...
        predLabelsText, ...
        procrustesMap(key), ...
        TP, ...
        FP, ...
        FN, ...
        TN, ...
        sensitivity, ...
        specificity ...
    };

end

resultsTable = cell2table(outRows, ...
    'VariableNames', { ...
        'ImageName', ...
        'PcellLabels', ...
        'LineagePositiveLabels', ...
        'CeSCALEPositiveLabels_AdjustedToOriginalIDs', ...
        'Procrustes', ...
        'TP', ...
        'FP', ...
        'FN', ...
        'TN', ...
        'Sensitivity', ...
        'Specificity' ...
    });

% Write as .xlsx file
writetable(resultsTable, ...
    'CeSCALE_Lineage_Classification_Results.xlsx', ...
    'FileType', 'spreadsheet');

%Function used for both P-Cell and Lineage CSVs
function labelMap = makeLabelMapFromCSV(inputFile)

    % Reads a CSV where:
    % column 1 = image/mask filename
    % column 2 = comma-separated label list, e.g. "12,13"
    %
    % Returns:
    % labelMap('405_0_11_FOV1') = [12 13]

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
    csvLabelText = strings(length(lines)-1, 1);

    % Parse CSV lines, skipping header
    for i = 2:length(lines)

        firstComma = strfind(lines{i}, ',');

        if isempty(firstComma)
            continue
        end

        csvImageNames(i-1) = strtrim(lines{i}(1:firstComma(1)-1));
        csvLabelText(i-1) = strtrim(lines{i}(firstComma(1)+1:end));

    end

    % Clean parsed values
    csvImageNames = regexprep(csvImageNames, "_cp_masks.*", "");
    csvLabelText = strrep(csvLabelText, '"', '');

    % Create dictionary: image name -> numeric label list
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