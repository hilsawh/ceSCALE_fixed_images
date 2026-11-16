%% CeSCALE confidence overlays for one embryo (1-P-cell datasets)
% Overlays CeSCALE cell identities and per-cell confidence scores on 3D
% median-filtered (3x3x3 voxels) maximum-intensity projections of the DAPI,
% P granule and ABa lineage-marker channels.
%
% Inputs (under dataRoot):
%   combined_per_cell_confidence.csv             per-cell CeSCALE confidence
%   Ncoords_1/parsed/Sample<N>FOV<M>coords.mat   crop coordinates (x1 x2 y1 y2)
%   Cropped_Embryo_Images/DAPI_0_<N>_FOV<M>.mat  cropped DAPI stack (I3)
%   hyb0_channels/PGranule_0_<N>.mat             P granule stack
%   hyb0_channels/Lineage_0_<N>.mat              ABa lineage-marker stack
% Plus, in the current folder:
%   viridis256.csv                               viridis (matplotlib v3.11.2)
%
% Outputs (in overlays/): .tif + .pdf per panel, colour bar .pdf
% Requires: MATLAB R2020a+ (exportgraphics), Image Processing Toolbox (medfilt3)

clearvars
close all

%% 0. Settings (edit this section only)
dataRoot  = '/path/to/data';   % folder with CeSCALE outputs + image stacks
runTag    = "dataset_2_R2";      % run prefix of sample_id in the confidence file
SampleNum = 15;
FOV       = 1;

pixelSize         = 0.11;       % um per pixel
scaleBarLength_um = 5;          % um
medK              = [3 3 3];    % medfilt3 neighbourhood (rows, cols, planes); must be odd
confRange         = [0 0.7];    % colour-bar range; scores >= 0.7 get the top colour
cmapFile          = 'viridis256.csv';

dotSize   = 40;          % confidence dots, DAPI panel
abaColor  = '#F2AF0D';   % ABa markers, lineage panel
zColor    = '#00FF00';   % P-cell marker, P granule panel
germColor = [1 1 1];     % P-cell marker, DAPI panel (white: germline cells are anchored and not scored)
noScoreColor = [0.5 0.5 0.5];   % grey: cells without a confidence score

%% 1. Setup
embryoid = ['DAPI_0_' num2str(SampleNum) '_FOV' num2str(FOV)];
confID   = runTag + "_" + string(embryoid);    % sample_id in the confidence file, e.g. dataset_2_R2_DAPI_0_15_FOV1

cmapConf = readmatrix(cmapFile);
assert(isequal(size(cmapConf), [256 3]), '%s should be 256 x 3', cmapFile)
% confidence -> RGB; scores outside confRange clip to the end colours
conf2rgb = @(c) cmapConf(round(min(max((c - confRange(1)) / diff(confRange), 0), 1) * 255) + 1, :);

tickVals = [0 0.25 0.5 0.75 1];
tickVals = unique([tickVals(tickVals >= confRange(1) & tickVals < confRange(2)), confRange(2)]);

save_path = 'overlays/';                 % relative to the current folder
if ~isfolder(save_path), mkdir(save_path), end
panelStem = [save_path char(runTag) '_Sample' num2str(SampleNum) 'FOV' num2str(FOV)];

%% 2. Load data and make median-filtered max projections
confFile = fullfile(dataRoot, 'combined_per_cell_confidence.csv');
confOpts = detectImportOptions(confFile);
confOpts = setvartype(confOpts, {'sample_id', 'imageID', 'Assigned_name'}, 'string');
conf = readtable(confFile, confOpts);
% MATLAB renames the "Relative_self_confidence(Procrustes_mod)" header, so find it by prefix
conf.confidence = conf{:, startsWith(conf.Properties.VariableNames, 'Relative_self_confidence')};
clear confOpts

crop = load(fullfile(dataRoot, 'Ncoords_1', 'parsed', sprintf('Sample%dFOV%dcoords.mat', SampleNum, FOV)));
x1 = crop.x1;  x2 = crop.x2;  y1 = crop.y1;  y2 = crop.y2;

D = load(fullfile(dataRoot, 'Cropped_Embryo_Images', [embryoid '.mat']));   % I3, already cropped
ImDAPI = max(medfilt3(D.I3, medK), [], 3);
clear D

S = load(fullfile(dataRoot, 'hyb0_channels', sprintf('PGranule_0_%d.mat', SampleNum)));
v = struct2cell(S); assert(isscalar(v), 'Variables in file: %s', strjoin(fieldnames(S), ', '))
Im561 = max(medfilt3(v{1}(y1:y2, x1:x2, :), medK), [], 3);

S = load(fullfile(dataRoot, 'hyb0_channels', sprintf('Lineage_0_%d.mat', SampleNum)));
v = struct2cell(S); assert(isscalar(v), 'Variables in file: %s', strjoin(fieldnames(S), ', '))
Im488 = max(medfilt3(v{1}(y1:y2, x1:x2, :), medK), [], 3);
clear S v

xlims = [1 size(ImDAPI, 2)];
ylims = [1 size(ImDAPI, 1)];

%% 3. Select this embryo's cells
embConf = conf(conf.sample_id == confID, :);
assert(height(embConf) > 0, 'No confidence scores for %s (embryo may not have been scored).', confID)

% germline cells (Z2/Z3 or a P cell) are anchored, not scored
shortNames = extractBefore(embConf.Assigned_name, min(4, strlength(embConf.Assigned_name) + 1));   % first 3 characters
isP    = startsWith(shortNames, "P");
isGerm = contains(shortNames, 'Z2') | contains(shortNames, 'Z3') | isP;
pXY    = [embConf.Initial_X(isP) embConf.Initial_Y(isP)];   % P-cell position(s)
pNames = cellstr(embConf.Assigned_name(isP));
embConf = embConf(~isGerm, :);

isABa = startsWith(embConf.Assigned_name, "ABa");

dotRGB = conf2rgb(embConf.confidence);
isNoScore = isnan(embConf.confidence);   % NaN would otherwise silently take the bottom colour
if any(isNoScore)
    warning('%d cell(s) in %s have no confidence score; shown in grey.', sum(isNoScore), confID)
    dotRGB(isNoScore, :) = repmat(noScoreColor, sum(isNoScore), 1);
end

Procrustes = embConf.min_Procrustes_mod_distance(1);   % one value per embryo; used in titles

%% 4. DAPI, nuclei coloured by CeSCALE confidence
figure('Renderer', 'painters', 'Position', [50 50 400 400])
imagesc(ImDAPI)
colormap gray
hold on
scatter(embConf.Initial_X, embConf.Initial_Y, dotSize, dotRGB, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.25)
if ~isempty(pXY), scatter(pXY(:,1), pXY(:,2), 25, germColor, 'd', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.25), end
title([embryoid ' DAPI pr=' num2str(Procrustes)])
axis equal
xlim(xlims)
ylim(ylims)
set(gca, 'FontSize', 7)

% scale bar (bottom-right), drawn in pixel units
barLength_px = scaleBarLength_um / pixelSize;
xEnd = xlims(2) - 0.05 * diff(xlims);
yBar = ylims(2) - 0.05 * diff(ylims);   
plot([xEnd - barLength_px, xEnd], [yBar yBar], 'w-', 'LineWidth', 3)

saveas(gcf, [panelStem '_DAPI_confidence.tif'])
exportgraphics(gcf, [panelStem '_DAPI_confidence.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');

%% 5. P granule stain, P cell(s)
figure('Renderer', 'painters', 'Position', [50 50 400 400])
imagesc(Im561)
colormap gray
hold on
scatter(pXY(:,1), pXY(:,2), 40, 'd', 'filled', 'MarkerFaceColor', zColor)
text(pXY(:,1), pXY(:,2) + 5, pNames, 'Color', 'white', 'FontSize', 7)
title([embryoid ' P granule pr=' num2str(Procrustes)])
axis equal
xlim(xlims)
ylim(ylims)
set(gca, 'FontSize', 7)
saveas(gcf, [panelStem '_PGranule.tif'])
exportgraphics(gcf, [panelStem '_PGranule.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');

%% 6. ABa lineage stain, CeSCALE ABa cells
figure('Renderer', 'painters', 'Position', [50 50 400 400])
imagesc(Im488)
colormap gray
hold on
scatter(embConf.Initial_X(isABa), embConf.Initial_Y(isABa), 25, 'filled', 'MarkerFaceColor', abaColor)
text(embConf.Initial_X(isABa), embConf.Initial_Y(isABa), cellstr(embConf.Assigned_name(isABa)), 'Color', 'white', 'FontSize', 7)
title([embryoid ' Lineage - ABa cells (n=' num2str(sum(isABa)) ')'])
axis equal
xlim(xlims)
ylim(ylims)
set(gca, 'FontSize', 7)
saveas(gcf, [panelStem '_Lineage_ABa.tif'])
exportgraphics(gcf, [panelStem '_Lineage_ABa.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');

%% 7. Colour bar for the confidence dots
figure('Renderer', 'painters', 'Position', [50 50 300 120])
axCb = axes('Visible', 'off');
colormap(axCb, cmapConf)
set(axCb, 'CLim', confRange)
cb = colorbar(axCb, 'southoutside');
labels = string(tickVals);
if confRange(2) < 1, labels(end) = "\geq" + confRange(2); end   % capped: top reads >= 0.7
cb.Ticks = tickVals;
cb.TickLabels = labels;
cb.Label.String = 'CeSCALE confidence';
cb.FontSize = 7;
exportgraphics(gcf, [save_path 'confidence_colorbar.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');