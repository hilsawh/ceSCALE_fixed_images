% create a signal file of transcripts in particular cells binned into
% specified genomic locations

clear all
close all

pat = 'D';
individual_cell_analysis = 0; %decide if looking at individual cell using strcmp (1) or entire lineage using 'startsWith' (0)

%% Import scRNAseq data from Cole et al eLife

opts = delimitedTextImportOptions("NumVariables", 122, "Encoding", "UTF-8");

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Var1", "Var2", "GeneStableID", "P0", "AB", "P1", "ABa", "ABp", "EMS", "P2", "ABax", "ABpx", "MS", "E", "C", "P38cellStage", "ABala", "ABalp", "ABara", "ABarp", "ABpxa", "ABpxp", "MSx", "Ea_15c", "Ep_15c", "Cx", "P315", "ABalaa", "ABalap", "ABalpa", "ABalpp", "ABaraa", "ABarap", "ABarpa", "ABarpp", "ABpxaa", "ABpxap", "ABpxpa", "ABpxpp", "MSxa", "MSxp", "Ea_28c", "Ep_28c", "Cxa", "Cxp", "D", "P428cellStage", "ABalaaa", "ABalaap", "ABalapa", "ABalapp", "ABalpaa", "ABalpap", "ABalppa", "ABalppp", "ABaraaa", "ABaraap", "ABarapa", "ABarapp", "ABarpaa", "ABarppa", "ABpxaaa", "ABpxaap", "ABpxapa", "ABpxapp", "ABpxpaa", "ABpxpap", "ABpxppa", "ABpxppp", "MSaaa", "MSaap", "MSapa", "MSapp", "MSpaa", "MSpap", "MSppa", "MSppp", "Eax", "Epx", "Caa_51c", "Cap_51c", "Cpa_51c", "Cpp_51c", "Dx", "P451cellStage", "ABalaaax", "ABalaapx", "ABalapax", "ABalappx", "ABalpaax", "ABalpapx", "ABalppax", "ABalpppx", "ABaraaax", "ABaraapx", "ABarapax", "ABarappx", "ABarpaxx", "ABarppax", "ABarpppx", "ABpxaaax", "ABpxaapx", "ABpxapax", "ABpxappx", "ABpxpaax", "ABpxpapx", "ABpxppax", "ABpxpppx", "MSxaaa", "MSxaap", "MSxapa", "MSxpaa", "MSxppa", "MSxppp", "Eala102", "Epla102", "Caaa", "Capa", "Cpaa", "Cppa", "Daa", "P4a102"];
opts.SelectedVariableNames = ["GeneStableID", "P0", "AB", "P1", "ABa", "ABp", "EMS", "P2", "ABax", "ABpx", "MS", "E", "C", "P38cellStage", "ABala", "ABalp", "ABara", "ABarp", "ABpxa", "ABpxp", "MSx", "Ea_15c", "Ep_15c", "Cx", "P315", "ABalaa", "ABalap", "ABalpa", "ABalpp", "ABaraa", "ABarap", "ABarpa", "ABarpp", "ABpxaa", "ABpxap", "ABpxpa", "ABpxpp", "MSxa", "MSxp", "Ea_28c", "Ep_28c", "Cxa", "Cxp", "D", "P428cellStage", "ABalaaa", "ABalaap", "ABalapa", "ABalapp", "ABalpaa", "ABalpap", "ABalppa", "ABalppp", "ABaraaa", "ABaraap", "ABarapa", "ABarapp", "ABarpaa", "ABarppa", "ABpxaaa", "ABpxaap", "ABpxapa", "ABpxapp", "ABpxpaa", "ABpxpap", "ABpxppa", "ABpxppp", "MSaaa", "MSaap", "MSapa", "MSapp", "MSpaa", "MSpap", "MSppa", "MSppp", "Eax", "Epx", "Caa_51c", "Cap_51c", "Cpa_51c", "Cpp_51c", "Dx", "P451cellStage", "ABalaaax", "ABalaapx", "ABalapax", "ABalappx", "ABalpaax", "ABalpapx", "ABalppax", "ABalpppx", "ABaraaax", "ABaraapx", "ABarapax", "ABarappx", "ABarpaxx", "ABarppax", "ABarpppx", "ABpxaaax", "ABpxaapx", "ABpxapax", "ABpxappx", "ABpxpaax", "ABpxpapx", "ABpxppax", "ABpxpppx", "MSxaaa", "MSxaap", "MSxapa", "MSxpaa", "MSxppa", "MSxppp", "Eala102", "Epla102", "Caaa", "Capa", "Cpaa", "Cppa", "Daa", "P4a102"];
opts.VariableTypes = ["string", "string", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["Var1", "Var2", "GeneStableID"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "Var2", "GeneStableID"], "EmptyFieldRule", "auto");

% Import the data
AverageTranscriptsPerCell = readtable("/scicore/home/mangos/GROUP/Kostas_Ahilya_shared_folder/Cole/AverageTranscriptsPerCell.csv", opts);

clear opts
%% Import gene locations from ParaSite Biomart using GeneStableID

opts = delimitedTextImportOptions("NumVariables", 5, "Encoding", "UTF-8");
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["GeneStableID", "GeneName", "ChromosomescaffoldName", "GeneStartbp", "GeneEndbp"];
opts.VariableTypes = ["string", "string", "categorical", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts = setvaropts(opts, ["GeneStableID", "GeneName"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["GeneStableID", "GeneName", "ChromosomescaffoldName"], "EmptyFieldRule", "auto");

TranscriptLocations = readtable("mart_export_GeneStableIDs.csv", opts);
clear opts
%% join tables by column name to connect transcript signal to gene location 
T1 = join(TranscriptLocations,AverageTranscriptsPerCell,'Keys','GeneStableID');

%% select columns of table based on name

%get all scRNA data column names with pat and find locations 
columns = T1.Properties.VariableNames(:,7:end); %exclude the gene name etc and P0 (maternal load) columns 1-6
if individual_cell_analysis ==0
    idx = find(startsWith(columns, pat));
    mkdir TXsignal_lineages_zygotic2
    save_dir = 'TXsignal_lineages_zygotic2/';
else
    idx = find(strcmp(columns, pat));
    mkdir TXsignal_cells
    save_dir = 'TXsignal_cells/';
end

%make table with tx signal from only cells containing pat
idx2 = idx +6;
T2 = T1(:, idx2);
names = T2.Properties.VariableNames;
T = T1(:, 1:5);

%% load zygotic gene list from clustering

opts = delimitedTextImportOptions("NumVariables", 1);
opts.DataLines = [1, Inf];
opts.Delimiter = ",";
opts.VariableNames = "WBGene00007063";
opts.VariableTypes = "string";
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts = setvaropts(opts, "WBGene00007063", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "WBGene00007063", "EmptyFieldRule", "auto");
zygotic_genes = readmatrix("non_maternal_expression_list.txt", opts);

clear opts
%% make chrI signal for each cell in pat and calculate variance

chromosome = 'I';
numregions = 14;

load(['TracingCoordschr' chromosome '.mat']) %chromosome bins

for k = 1:length(idx)
    %find genes in chromosome of interest
    LocList = T.GeneStartbp(T.ChromosomescaffoldName == chromosome,:);
    GeneList = T.GeneStableID(T.ChromosomescaffoldName == chromosome);
    GeneNames = T.GeneName(T.ChromosomescaffoldName == chromosome);
    CellSignalList = T2.(names{1,k})(T.ChromosomescaffoldName == chromosome);

    %assign each gene to closest tracing bin
    TXperregionI = zeros(numregions,1);  
    TracingCoordschrI.probe_centroid = (TracingCoordschrI.probestart + TracingCoordschrI.probestop)/2;
    BinAssignmentList = zeros(length(GeneList),1);

    for i = 1:length(GeneList)
        [~,closestRegion] = min(abs(LocList(i,1)-TracingCoordschrI.probe_centroid));%find closest probe_centroid to gene
        BinAssignmentList(i,1) = closestRegion;
    end

    %count signal intensity per bin for zygotic genes only
    for i = 1:numregions
        RegionList = [];
        NameList = [];
        for j = 1:length(GeneList)
            if BinAssignmentList(j,1)==i && ismember(GeneList(j,1), zygotic_genes)
                TXperregionI(i,1)= TXperregionI(i,1) + CellSignalList(j,1);
                RegionList = [RegionList TXperregionI(i,1)];
                NameList = [NameList GeneNames(j,1)];
            end            
        end
        RegionListI{i,k} = mean(RegionList);
        NameListI{i,k} = NameList;
    end
    
    % figure('Renderer', 'painters', 'Position', [450 650 300 200])
    % bar(TXperregionI)
    % xlabel(['chr' chromosome ' location'])
    % ylabel([ names{1,k} ' cell expression ' pat])
    % ylim([0 20000])
    % set(gca,'fontsize', 7);
    % set(gca, 'FontName', 'Arial');    
end
disp(['number of zygotic genes on chrI = ' num2str(sum(cellfun(@numel, NameListI)))])

StdRegionListI = [];
MeanRegionListI = [];
for i = 1:numregions
    StdRegionListI(i,1) = std(cell2mat(RegionListI(i,:)));
    MeanRegionListI(i,1) = mean(cell2mat(RegionListI(i,:)));
end

% figure('Renderer', 'painters', 'Position', [450 650 300 200])
% h = heatmap(cell2mat(RegionListI));
% h.XDisplayLabels = names;
% h.ColorLimits = [-max(max(abs(cell2mat(RegionListI)))) max(max((abs(cell2mat(RegionListI)))))];
% colormap redbluecmap
% ylabel(['chr' chromosome ' location'])
% xlabel([ 'individual cell expression ' pat])
% set(gca,'fontsize', 7);
% set(gca, 'FontName', 'Arial');

% figure('Renderer', 'painters', 'Position', [450 650 300 200])
% bar(MeanRegionListI, 'FaceColor', 'k')
% hold on
% er = errorbar(1:numregions, MeanRegionListI,StdRegionListI, StdRegionListI);
% er.LineStyle = 'none';
% er.Color = [0 0 0];
% xlabel(['chr' chromosome ' location'])
% ylabel(['average cell expression ' pat])
% %ylim([0 20000])
% set(gca,'fontsize', 7);
% set(gca, 'FontName', 'Arial');

clear chromosome LocList GeneNames CellSignalList GeneList BinAssignmentList

%% make chrII siganl for each cell in pat and calculate variance
chromosome = 'II';
load(['TracingCoordschr' chromosome '.mat']) %chromosome bins

for k = 1:length(idx)
    %find genes in chromosome of interest
    LocList = T.GeneStartbp(T.ChromosomescaffoldName == chromosome,:);
    GeneList = T.GeneStableID(T.ChromosomescaffoldName == chromosome);
    GeneNames = T.GeneName(T.ChromosomescaffoldName == chromosome);
    CellSignalList = T2.(names{1,k})(T.ChromosomescaffoldName == chromosome);

    %assign each gene to closest tracing bin
    TXperregionII = zeros(numregions,1);  
    TracingCoordschrII.probe_centroid = (TracingCoordschrII.probestart + TracingCoordschrII.probestop)/2;
    BinAssignmentList = zeros(length(GeneList),1);

    for i = 1:length(GeneList)
        [~,closestRegion] = min(abs(LocList(i,1)-TracingCoordschrII.probe_centroid));%find closest probe_centroid to gene
        BinAssignmentList(i,1) = closestRegion;
    end

    %count signal intensity per bin
    for i = 1:numregions
        RegionList = [];
        NameList = [];
        for j = 1:length(GeneList)
            if BinAssignmentList(j,1)==i && ismember(GeneList(j,1), zygotic_genes)
                TXperregionII(i,1)= TXperregionII(i,1) + CellSignalList(j,1);
                RegionList = [RegionList TXperregionII(i,1)];
                NameList = [NameList GeneNames(j,1)];
            end            
        end
        RegionListII{i,k} = mean(RegionList);
        NameListII{i,k} = NameList;
    end

    % figure('Renderer', 'painters', 'Position', [450 650 300 200])
    % bar(TXperregionI)
    % xlabel(['chr' chromosome ' location'])
    % ylabel([ names{1,k} ' cell expression ' pat])
    % ylim([0 20000])
    % set(gca,'fontsize', 7);
    % set(gca, 'FontName', 'Arial');    
end
disp(['number of zygotic genes on chrII = ' num2str(sum(cellfun(@numel, NameListII)))])
StdRegionListII = [];
MeanRegionListII = [];
for i = 1:numregions
    StdRegionListII(i,1) = std(cell2mat(RegionListII(i,:)));
    MeanRegionListII(i,1) = mean(cell2mat(RegionListII(i,:)));
end

% figure('Renderer', 'painters', 'Position', [450 650 300 200])
% h = heatmap(cell2mat(RegionListII));
% colormap redbluecmap
% h.XDisplayLabels = names;
% h.ColorLimits = [-max(max(abs(cell2mat(RegionListI)))) max(max((abs(cell2mat(RegionListI)))))];
% ylabel(['chr' chromosome ' location'])
% xlabel([ 'individual cell expression ' pat])
% set(gca,'fontsize', 7);
% set(gca, 'FontName', 'Arial');

% figure('Renderer', 'painters', 'Position', [450 650 300 200])
% bar(MeanRegionListII, 'FaceColor', 'k')
% hold on
% er = errorbar(1:numregions, MeanRegionListII,StdRegionListII, StdRegionListII);
% er.LineStyle = 'none';
% er.Color = [0 0 0];
% xlabel(['chr' chromosome ' location'])
% ylabel(['average cell expression ' pat])
% %ylim([0 20000])
% set(gca,'fontsize', 7);
% set(gca, 'FontName', 'Arial');


%% plot merged signal

RegionListTotal = vertcat(RegionListI, RegionListII);
figure('Renderer', 'painters', 'Position', [650 150 300 400])
h = heatmap(cell2mat(RegionListTotal));
colormap copper
h.XDisplayLabels = names;
%h.ColorLimits = [-max(max(abs(cell2mat(RegionListTotal)))) max(max((abs(cell2mat(RegionListTotal)))))];
h.ColorLimits = [0 10000];
ylabel('region')
xlabel([ 'individual cell expression ' pat])
set(gca,'fontsize', 7);
%set(gca, 'FontName', 'Arial');

exportgraphics(gcf, [save_dir '/MeanTXsignalHeatmap_ ' pat '.pdf'],'ContentType','vector','BackgroundColor','none');


MeanRegionListTotal = vertcat(MeanRegionListI, MeanRegionListII);
StdRegionListTotal = vertcat(StdRegionListI, StdRegionListII);
NameListTotal = vertcat(NameListI, NameListII);

figure('Renderer', 'painters', 'Position', [450 650 300 200])
bar(MeanRegionListTotal, 'FaceColor', 'k')
hold on
er = errorbar(1:28, MeanRegionListTotal,StdRegionListTotal, StdRegionListTotal);
er.LineStyle = 'none';
er.Color = [0 0 0];
xlabel('region')
ylabel(['average expression ' pat ' lineage'])
%ylim([0 20000])
set(gca,'fontsize', 7);
%set(gca, 'FontName', 'Arial');

exportgraphics(gcf, [save_dir '/MeanTXsignalBar_ ' pat '.pdf'],'ContentType','vector','BackgroundColor','none');

logMeanRegionListTotal = log10(MeanRegionListTotal);

figure('Renderer', 'painters', 'Position', [650 150 300 400])
tiledlayout(1,2)
ax1=nexttile;
h = heatmap(MeanRegionListTotal);
colormap copper
%h.XDisplayLabels = names;
%h.ColorLimits = [-max(max(abs(cell2mat(RegionListTotal)))) max(max((abs(cell2mat(RegionListTotal)))))];
h.ColorLimits = [0 10000];
ylabel('region')
xlabel([ 'expression ' pat])
set(gca,'fontsize', 7);

ax2=nexttile;
h = heatmap(logMeanRegionListTotal);
colormap copper
%h.XDisplayLabels = names;
%h.ColorLimits = [-max(max(abs(cell2mat(RegionListTotal)))) max(max((abs(cell2mat(RegionListTotal)))))];
%h.ColorLimits = [0 10000];
ylabel('region')
xlabel([ 'log10 expression ' pat])
set(gca,'fontsize', 7);

exportgraphics(gcf, [save_dir '/RawvsLogMeanExpression_ ' pat '.pdf'],'ContentType','vector','BackgroundColor','none');


%% save the output for MetaLoci spatial autocorrelation

csvwrite([save_dir 'MeanTxsignal_' pat '.csv'], MeanRegionListTotal)
csvwrite([save_dir 'log10MeanTxsignal_' pat '.csv'], logMeanRegionListTotal)
csvwrite([save_dir 'StdTxsignal_' pat '.csv'], StdRegionListTotal)
save([save_dir 'GenesTxsignal_' pat '.mat'], 'NameListTotal')


