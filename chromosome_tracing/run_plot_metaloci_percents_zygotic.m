clear all
close all

type = 'lh';

lineage_list = {'AB', 'MS', 'E', 'C', 'D'};

mkdir('figures_zygotic')
save_dir = 'figures_zygotic';

for i = 1:length(lineage_list)
    lineage = lineage_list{1,i};
    %% import metaloci percentages data
    opts = delimitedTextImportOptions("NumVariables", 31);

    opts.DataLines = [2, Inf];
    opts.Delimiter = ",";

    opts.VariableNames = ["VarName1", "VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11", "VarName12", "VarName13", "VarName14", "VarName15", "VarName16", "VarName17", "VarName18", "VarName19", "VarName20", "VarName21", "VarName22", "VarName23", "VarName24", "VarName25", "VarName26", "VarName27", "VarName28", "VarName29", "traces", "cell_id"];
    opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string"];

    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

    opts = setvaropts(opts, "cell_id", "WhitespaceRule", "preserve");
    opts = setvaropts(opts, "cell_id", "EmptyFieldRule", "auto");

    data = readtable(['Percentages_zygotic/' type '_data_' lineage '_zygotic.csv'], opts);

    clear opts
    %% plot heatmap and line graph of mean + sd

    names = data.cell_id(1:end);
    percents = data{1:end,2:29};
    mean_percents = mean(percents);
    std_percents = std(percents);

    figure('Renderer', 'painters', 'Position', [50 50 600 800])
    subplot('Position', [0.15 0.8 0.5 0.1])
    hold on
    stdshade(percents, 0.9, [0.9 0.9 0.9])
    plot(mean_percents,'Color', 'k', 'LineWidth',1,'LineStyle','--')
    xlim([0.5 28.5])
    %ylim([0 0.3])
    set(gca, 'XTickLabel', [])
    set(gca,'fontsize', 7);
    %set(gca, 'FontName', 'Arial');
    ylabel('percentage')
    title([lineage ' lineage ' type ' metaloci'])

    subplot('Position', [0.15 0.35 0.5 0.4])
    h = heatmap(percents);
    colormap(flipud(bone))
    clim([0 0.4])
    h.YData = [names];
    grid off
    s = struct(h);
    s.XAxis.TickLabelRotation = 45;
    set(gca,'fontsize', 7);
    %set(gca, 'FontName', 'Arial');
    xlabel('region')
    ylabel('cell ID')

    exportgraphics(gcf, [save_dir '/metaloci_heatmap_' lineage '_' type '.pdf'],'ContentType','vector','BackgroundColor','none');
    csvwrite([save_dir '/mean_percents_' lineage '_' type '.csv'], mean_percents);

    %% plot clustergram https://www.mathworks.com/help/bioinfo/ug/working-with-the-clustergram-function.html

    cg_s = clustergram(percents, 'RowLabels', names,...
        'Cluster', 'Column',...
        'ImputeFun', @knnimpute)

    cg_s.ColumnPDist = 'euclidean';
    cg_s.Colormap = redbluecmap;
    set(0,'ShowHiddenHandles','on')
    % Get all handles from root
    allhnds = get(0,'Children');
    % Find hearmap axis and change the font size
    h = findall(allhnds, 'Tag', 'HeatMapAxes');
    set(h, 'FontSize', 7)

    set(gcf, 'Position', [550 250 600 500])

    exportgraphics(gcf, [save_dir '/metaloci_clustergram_' lineage '_' type '.pdf'],'ContentType','vector','BackgroundColor','none');

end

















