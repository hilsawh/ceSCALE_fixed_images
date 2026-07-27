% create a signal file of transcripts in particular cells binned into
% specified genomic locations

clear all
close all

%% load cell names
opts = delimitedTextImportOptions("NumVariables", 2, "Encoding", "UTF-8");
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["Var1", "ColeCellTranscriptSignalCellName"];
opts.SelectedVariableNames = "ColeCellTranscriptSignalCellName";
opts.VariableTypes = ["string", "string"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts = setvaropts(opts, ["Var1", "ColeCellTranscriptSignalCellName"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "ColeCellTranscriptSignalCellName"], "EmptyFieldRule", "auto");

MetaLociCellsV2 = readmatrix("MetaLociCellsV2.csv", opts);

clear opts

%% Import scRNAseq data from Cole et al eLife

opts = delimitedTextImportOptions("NumVariables", 120, "Encoding", "UTF-8");

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types %%%%% correct ABxpap to ABpxpap, typo in Cole paper
opts.VariableNames = ["GeneStableID", "P0", "AB", "P1", "ABa", "ABp", "EMS", "P2", "ABax", "ABpx", "MS", "E", "C", "P3_8c", "ABala", "ABalp", "ABara", "ABarp", "ABpxa", "ABpxp", "MSx", "Ea_15c", "Ep_15c", "Cx", "P3_15c", "ABalaa", "ABalap", "ABalpa", "ABalpp", "ABaraa", "ABarap", "ABarpa", "ABarpp", "ABpxaa", "ABpxap", "ABpxpa", "ABpxpp", "MSxa", "MSxp", "Ea_28c", "Ep_28c", "Cxa", "Cxp", "D", "P4_28c", "ABalaaa", "ABalaap", "ABalapa", "ABalapp", "ABalpaa", "ABalpap", "ABalppa", "ABalppp", "ABaraaa", "ABaraap", "ABarapa", "ABarapp", "ABarpaa", "ABarppa", "ABpxaaa", "ABpxaap", "ABpxapa", "ABpxapp", "ABpxpaa", "ABpxpap", "ABxrppa", "ABxrppp", "MSaaa", "MSaap", "MSapa", "MSapp", "MSpaa", "MSpap", "MSppa", "MSppp", "Eax", "Epx", "Caa_51c", "Cap_51c", "Cpa_51c", "Cpp_51c", "Dx", "P4_51c", "ABalaaax", "ABalaapx", "ABalapax", "ABalappx", "ABalpaax", "ABalpapx", "ABalppax", "ABalpppx", "ABaraaax", "ABaraapx", "ABarapax", "ABarappx", "ABarpaxx", "ABarppax", "ABarpppx", "ABpxaaax", "ABpxaapx", "ABpxapax", "ABpxappx", "ABpxpaax", "ABpxpapx", "ABpxppax", "ABpxpppx", "MSxaaa", "MSxaap", "MSxapa", "MSxpaa", "MSxppa", "MSxppp", "Eala_102c", "Epla_102c", "Caaa", "Capa", "Cpaa", "Cppa", "Daa", "P4_102c"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "GeneStableID", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "GeneStableID", "EmptyFieldRule", "auto");

% Import the data
AverageTranscriptsPerCell = readtable("AverageTranscriptsPerCell.csv", opts);

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
mkdir TXsignal_cells_zygotic
save_dir = 'TXsignal_cells_zygotic/';

for q = 1:length(MetaLociCellsV2)
    pat = MetaLociCellsV2(q, 1);
    %get all scRNA data column names with pat and find locations
    columns = T1.Properties.VariableNames(:,7:end); %exclude the gene name etc and P0 (maternal load) columns 1-6
    idx = find(strcmp(columns, pat));

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

        %count signal intensity per bin for zygotic genes
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
    end
    disp(['number of zygotic genes on chrI = ' num2str(sum(cellfun(@numel, NameListI)))])
    
    StdRegionListI = [];
    MeanRegionListI = [];
    for i = 1:numregions
        StdRegionListI(i,1) = std(cell2mat(RegionListI(i,:)));
        MeanRegionListI(i,1) = mean(cell2mat(RegionListI(i,:)));
    end

clear chromosome LocList GeneNames CellSignalList GeneList BinAssignmentList

    %% make chrII signal for each cell in pat and calculate variance
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
        
    end
    disp(['number of zygotic genes on chrII = ' num2str(sum(cellfun(@numel, NameListII)))])

    StdRegionListII = [];
    MeanRegionListII = [];
    for i = 1:numregions
        StdRegionListII(i,1) = std(cell2mat(RegionListII(i,:)));
        MeanRegionListII(i,1) = mean(cell2mat(RegionListII(i,:)));
    end

   
    %% plot merged signal

    RegionListTotal = vertcat(RegionListI, RegionListII);
    figure('Renderer', 'painters', 'Position', [650 150 300 400])
    h = heatmap(cell2mat(RegionListTotal));
    colormap redbluecmap
    h.XDisplayLabels = names;
    h.ColorLimits = [-max(max(abs(cell2mat(RegionListTotal)))) max(max((abs(cell2mat(RegionListTotal)))))];
    ylabel(['chr location'])
    xlabel([ 'individual cell expression ' pat])
    set(gca,'fontsize', 7);
    %set(gca, 'FontName', 'Arial');

    exportgraphics(gcf, strcat(save_dir, "MeanTXsignalHeatmap_",  pat, ".pdf"),'ContentType','vector','BackgroundColor','none');
  

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
    ylabel(['average expression ' pat])
    %ylim([0 20000])
    set(gca,'fontsize', 7);
    %set(gca, 'FontName', 'Arial');

    exportgraphics(gcf, strcat(save_dir, 'MeanTXsignalBar_', pat, '.pdf'),'ContentType','vector','BackgroundColor','none');

    %% save the output

    csvwrite(strcat(save_dir, 'MeanTxsignal_', pat, '.csv'), MeanRegionListTotal)
    csvwrite(strcat(save_dir, 'StdTxsignal_', pat, '.csv'), StdRegionListTotal)
    save(strcat(save_dir, 'GenesTxsignal_', pat, '.mat'), 'NameListTotal')

end