%% CeSCALE confidence overlays for one embryo (2-P-cell datasets)
% Overlays CeSCALE cell identities and per-cell confidence scores on
% maximum-intensity projections of the DAPI, P granule and E-cell
% lineage-marker channels. The lineage panel projects only the planes
% eSlice +/- eHalfWidth. No filtering is applied here; the stacks were
% denoised upstream with Denoise.ai in NIS-Elements.
%
% Inputs:
%   under dataRoot:
%     combined_per_cell_confidence.csv           (per-cell CeSCALE confidence) 
%   under imageRoot:
%     Ncoords_1/DAPI_0_<N>_FOV<M>.mat            Z2/Z3 centroids (z2center, z3center)
%     Ncoords_1/parsed/Sample<N>FOV<M>coords.mat   crop coordinates (x1 x2 y1 y2)
%     hyb0_channels/DAPI_0_<N>.mat, hyb0_channels/Lineage_0_<N>.mat, hyb0_channels/PGranule_0_<N>.mat   image stacks
%   in the current folder:
%     viridis256.csv                                 exact viridis (matplotlib v3.11.2, _cm_listed.py), 256 x 3
%
% Outputs (in overlays/): .tif + .pdf per panel, colour bar .pdf
% Requires: MATLAB R2020a+ (exportgraphics)

clearvars
close all

%% 0. Settings (edit this section only)
dataRoot     = '/path/to/data';   % folder with the confidence CSV
imageRoot    = '.';    % folder with Ncoords_1/ and hyb0_channels/ ('.' = current folder)
confFileName = 'combined_per_cell_confidence.csv';
runTag       = "R2";   % run prefix of sample_id in the confidence file
namePrefix   = '405_0_';   % start of the embryo name in sample_id ('405_0_' or 'DAPI_0_')
SampleNum    = 6;
FOV          = 1;

eSlice     = 125;   % centre z-slice of lineage channel, chosen for this embryo
eHalfWidth = 10;    % eSlice +/- eHalfWidth planes

pixelSize         = 0.11;       % um per pixel
scaleBarLength_um = 5;          % um
confRange         = [0 0.7];    % colour-bar range; scores >= 0.7 get the top colour
cmapFile          = 'viridis256.csv';

dotSize      = 25;          % confidence dots, DAPI panel
eColor       = '#F2AF0D';   % E-cell markers, Lineage panel
zColor       = '#00FF00';   % Z2/Z3 markers, PGranule panel
germColor    = [1 1 1];     % Z2/Z3/P markers, DAPI panel (white: germline cells are anchored, not scored)
noScoreColor = [0.5 0.5 0.5];   % grey: cells without a confidence score

%% 1. Setup
embryoid = [namePrefix num2str(SampleNum) '_FOV' num2str(FOV)];
confID   = runTag + "_merged" + string(embryoid);   % sample_id in the confidence file, e.g. R2_merged405_0_6_FOV1

cmapConf = readmatrix(cmapFile);
assert(isequal(size(cmapConf), [256 3]), '%s should be 256 x 3', cmapFile)
% confidence -> RGB; scores outside confRange clip to the end colours
conf2rgb = @(c) cmapConf(round(min(max((c - confRange(1)) / diff(confRange), 0), 1) * 255) + 1, :);

tickVals = [0 0.25 0.5 0.75 1];
tickVals = unique([tickVals(tickVals >= confRange(1) & tickVals < confRange(2)), confRange(2)]);

save_path = 'overlays/';                 % relative to the current folder
if ~isfolder(save_path), mkdir(save_path), end
panelStem = [save_path char(runTag) '_Sample' num2str(SampleNum) 'FOV' num2str(FOV)];

%% 2. Load data and make max projections
confFile = fullfile(dataRoot, confFileName);
confOpts = detectImportOptions(confFile);
confOpts = setvartype(confOpts, {'sample_id', 'imageID', 'Assigned_name'}, 'string');
conf = readtable(confFile, confOpts);
% MATLAB renames the "Relative_self_confidence(Procrustes_mod)" header, so find it by prefix
conf.confidence = conf{:, startsWith(conf.Properties.VariableNames, 'Relative_self_confidence')};
clear confOpts

% crop coordinates and Z2/Z3 centres
G = load(fullfile(imageRoot, 'Ncoords_1', ['DAPI_0_' num2str(SampleNum) '_FOV' num2str(FOV) '.mat']), ...
         'z2center', 'z3center');
C = load(fullfile(imageRoot, 'Ncoords_1', 'parsed', sprintf('Sample%dFOV%dcoords.mat', SampleNum, FOV)));
x1 = C.x1;  x2 = C.x2;  y1 = C.y1;
if isfield(C, 'y2'), y2 = C.y2; else, y2 = C.y3; end   % new files: y2; old files: y3
z2center = G.z2center;  z3center = G.z3center;
clear G C

chDir = fullfile(imageRoot, 'hyb0_channels');
S = load(fullfile(chDir, sprintf('DAPI_0_%d.mat', SampleNum)), 'ImageStackDAPI');
ImDAPI = max(S.ImageStackDAPI(y1:y2, x1:x2, :), [], 3);

S = load(fullfile(chDir, sprintf('PGranule_0_%d.mat', SampleNum)), 'ImageStackPGranule');   % P granules
Im561 = max(S.ImageStackPGranule(y1:y2, x1:x2, :), [], 3);

S = load(fullfile(chDir, sprintf('Lineage_0_%d.mat', SampleNum)), 'ImageStackLineage');    % E cell stain
nZ488  = size(S.ImageStackLineage, 3);
eRange = max(1, eSlice - eHalfWidth) : min(nZ488, eSlice + eHalfWidth);
Im488  = max(S.ImageStackLineage(y1:y2, x1:x2, eRange), [], 3);
clear S

xlims = [1 size(ImDAPI, 2)];
ylims = [1 size(ImDAPI, 1)];

%% 3. Select this embryo's cells
embConf = conf(conf.sample_id == confID, :);
assert(height(embConf) > 0, 'No confidence scores for %s (embryo may not have been scored).', confID)

% germline cells (Z2/Z3 or a P cell) are anchored, not scored: split them off
shortNames = extractBefore(embConf.Assigned_name, min(4, strlength(embConf.Assigned_name) + 1));   % first 3 characters
isP    = startsWith(shortNames, "P");
isGerm = contains(shortNames, 'Z2') | contains(shortNames, 'Z3') | isP;
pXY    = [embConf.Initial_X(isP) embConf.Initial_Y(isP)];   % P-cell position(s)
embConf = embConf(~isGerm, :);

isE = startsWith(embConf.Assigned_name, "E");

dotRGB = conf2rgb(embConf.confidence);
isNoScore = isnan(embConf.confidence);   % NaN would otherwise silently take the bottom colour
if any(isNoScore)
    warning('%d cell(s) in %s have no confidence score; shown in grey.', sum(isNoScore), confID)
    dotRGB(isNoScore, :) = repmat(noScoreColor, sum(isNoScore), 1);
end

Procrustes = embConf.min_Procrustes_mod_distance(1);   % one value per embryo; used in titles

%% 4. Panel: DAPI, nuclei coloured by CeSCALE confidence
figure('Renderer', 'painters', 'Position', [50 50 400 400])
imagesc(ImDAPI)
colormap gray
hold on
scatter(embConf.Initial_X, embConf.Initial_Y, dotSize, dotRGB, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.25)
scatter(z2center(1), z2center(2), 25, germColor, 'd', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.25)
scatter(z3center(1), z3center(2), 25, germColor, 'd', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.25)
text(z2center(1), z2center(2) + 5, 'Z2', 'Color', 'black', 'FontSize', 4)
text(z3center(1), z3center(2) + 5, 'Z3', 'Color', 'black', 'FontSize', 4)
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

%% 5. Panel: P granules, Z2/Z3
figure('Renderer', 'painters', 'Position', [50 50 400 400])
imagesc(Im561)
colormap gray
hold on
scatter(z2center(1), z2center(2), 25, 'd', 'filled', 'MarkerFaceColor', zColor)
scatter(z3center(1), z3center(2), 25, 'd', 'filled', 'MarkerFaceColor', zColor)
text(z2center(1), z2center(2) + 5, 'Z2/3', 'Color', 'white', 'FontSize', 4)
text(z3center(1), z3center(2) + 5, 'Z2/3', 'Color', 'white', 'FontSize', 4)
title([embryoid ' PGranules pr=' num2str(Procrustes)])
axis equal
xlim(xlims)
ylim(ylims)
set(gca, 'FontSize', 7)
saveas(gcf, [panelStem '_Pgranules.tif'])
exportgraphics(gcf, [panelStem '_Pgranules.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');

%% 6. Panel: Lineage E-cell stain, CeSCALE E cells
figure('Renderer', 'painters', 'Position', [50 50 400 400])
imagesc(Im488)
colormap gray
hold on
scatter(embConf.Initial_X(isE), embConf.Initial_Y(isE), 25, 'filled', 'MarkerFaceColor', eColor)
text(embConf.Initial_X(isE), embConf.Initial_Y(isE), cellstr(embConf.Assigned_name(isE)), 'Color', 'white', 'FontSize', 7)
title([embryoid ' Lineage - E cells - slices ' num2str(eRange(1)) '-' num2str(eRange(end)) '/' num2str(nZ488)])
axis equal
xlim(xlims)
ylim(ylims)
set(gca, 'FontSize', 7)
saveas(gcf, [panelStem '_488_Ecells.tif'])
exportgraphics(gcf, [panelStem '_488_Ecells.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');

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