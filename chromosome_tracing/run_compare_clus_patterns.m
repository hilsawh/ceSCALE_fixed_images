clear all
close all

NumClusters = 11;
res = 'sub0res06_sub1res08';
%% set patterns to be counted 
TargetAge = 0; %currently not ensuring trace comes from specific age, but from identities found in this age (can also be found before and after)

%save_dir = res;

%pats = {'ABal' 'ABar' 'ABpl' 'ABpr' 'MS' 'E' 'C' 'D' 'P'};
%general lineages --- use startsWith!!

%pats = {'ABaraa', 'ABarap','ABarpa','ABarpp', 'ABalaa', 'ABalap', 'ABalpa','ABalpp','ABpraa','ABprpa','ABprap','ABprpp','ABplaa','ABplap','ABplpa','ABplpp', 'MSa', 'MSp','Ea', 'Ep', 'Ca', 'Cp', 'D','P4'}; 
%24cell

%pats = {'ABara' 'ABarp' 'ABala' 'ABalp' 'ABpra' 'ABprp' 'ABpla' 'ABplp' 'MSa' 'MSp' 'Ea' 'Ep' 'Ca' 'Cp' 'D' 'P4'}; 
%16cell

%pats = {'ABara' 'ABarp' 'ABala' 'ABalp' 'ABpra' 'ABprp' 'ABpla' 'ABplp' 'MSa' 'MSp' 'Ea' 'Ep' 'Ca' 'Cp' 'P3'}; 
%15cell

%pats = {'ABar' 'ABal' 'ABpr' 'ABpl' 'MS' 'E' 'C' 'P3'}; 
%8cell

%pats = {'ABa' 'ABp' 'MS' 'E' 'C' 'P3'}; 
%6cell

%pats = {'ABa' 'ABp' 'EMS', 'P2'}; 
%4cell

save_dir = 'AB';

%pats = {'C' 'Ca' 'Cp' 'Caa' 'Cap' 'Cpa' 'Cpp'}; %'Caaa' 'Caap' 'Capa' 'Capp' 'Cpaa' 'Cpap' 'Cppa' 'Cppp'}; %C lineage
%pats = {'D' 'Da' 'Dp'}; %D lineage
%pats = {'P2' 'P3' 'P4'}; %P lineage
%pats = {'EMS' 'E' 'Ea' 'Ep' 'Eal' 'Ear' 'Epl' 'Epr'}; %E lineage
%pats = {'EMS' 'MS' 'MSa' 'MSp' 'MSaa' 'MSap' 'MSpa' 'MSpp' 'MSaaa' 'MSaap' 'MSapa' 'MSapp' 'MSpaa' 'MSpap' 'MSppa' 'MSppp' }; %MS lineage
%%%%pats = {'ABa', 'ABp', 'ABar', 'ABpr', 'ABal', 'ABpl','ABara', 'ABpra', 'ABala', 'ABpla', 'ABarp', 'ABprp', 'ABalp', 'ABplp','ABaraa', 'ABpraa', 'ABalaa', 'ABplaa', 'ABarap', 'ABprap', 'ABalap', 'ABplap','ABarpa', 'ABprpa', 'ABalpa', 'ABplpa','ABarpp', 'ABprpp', 'ABalpp', 'ABplpp'}; %AB lineage
%% load cell ids by lineage from file

%AB pats
% 
opts = delimitedTextImportOptions("NumVariables", 1, "Encoding", "UTF-8");
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
opts.VariableNames = "AB";
opts.VariableTypes = "string";
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts = setvaropts(opts, "AB", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "AB", "EmptyFieldRule", "auto");
cellids = readtable("AB_cellids.csv", opts);
clear opts

pats = cellids.AB(1:end)';

%% counting by pat per cluster and save patlog

patlog = cell2table(cell(NumClusters,length(pats)),'VariableNames', pats);

for q = 1:NumClusters

    load([res '/ChrClus' num2str(q) '.mat'])

    for curr_pat = pats
        n = 0;
        for i = 1:length(ChrClus)
            if strcmp(ChrClus(i).CellID, curr_pat) %for exact pat
            %if startsWith(ChrClus(i).CellID, curr_pat) %for all in lineage after pat
                n = n+1;
            end
        end
        
        patlog.(char(curr_pat)){q} = n;
    end
end

save([save_dir '/proportions_' num2str(TargetAge) '-cell.mat'], 'patlog')

patlog1 = table2array(patlog(:,1:length(pats)));

%% normalize

%normalize to total num of all patterns in each CLUSTER;
patlognorm = zeros(NumClusters,length(pats));
for b = 1:NumClusters
    patlognorm(b,:) = cell2mat(patlog1(b,:))./sum(cell2mat(patlog1(b,:)));
end

n_clus = zeros (NumClusters, 1);
for b = 1:NumClusters
    n_clus(b,1) = sum(cell2mat(patlog1(b,:)));
end

%normalize to total num of all patterns in each CELL ID;
patlognorm2 = zeros(NumClusters,length(pats));
for b = 1:length(pats)
    patlognorm2(:,b) = cell2mat(patlog1(:,b))./sum(cell2mat(patlog1(:,b)));
end

n_cells = zeros (length(pats),1);
for b = 1:length(pats)
    n_cells(b,1) = sum(cell2mat(patlog1(:,b)));
end

%double normalize by cluster and cell id
patlogdbnorm = zeros(NumClusters, length(pats));
for b = 1:NumClusters
    patlogdbnorm(b,:) = cell2mat(patlog1(b,:))./sum(cell2mat(patlog1(b,:)));
end
for b = 1:length(pats)
    patlogdbnorm(:,b) = patlogdbnorm(:,b)./sum(patlogdbnorm(:,b));
end

%% plot figs and save

figure('Renderer', 'painters', 'Position', [50 50 400 300])
subplot('Position', [0.1 0.3 0.1 0.4])
scatter(ones(length(n_clus)), [1:length(n_clus)], flip(n_clus), 'k', 'filled')
xlim([0.9 1.1])
ylim([0.5 length(n_clus)+0.5])
axis off
subplot('Position', [0.25 0.3 0.5 0.4])
h = heatmap(patlognorm, 'MissingDataColor',[0.8 0.8 0.8]);
colormap(flipud(bone))
title(['normalized CLUSTER lineage proportions ' num2str(TargetAge) '-cell']);
caxis([0 0.3])
h.XData = [pats];
h.YData = [1:NumClusters];
grid off
s = struct(h);
s.XAxis.TickLabelRotation = 45;
set(gca,'fontsize', 7);
xlabel('cell ID')
ylabel('cluster ID')
%print([save_dir '/cellIDsHeatmapNorm ' num2str(TargetAge) ' cell'],'-dsvg');
exportgraphics(gcf, [save_dir '/cellIDsHeatmapNorm1 ' num2str(TargetAge) '-cell.pdf'],'ContentType','vector','BackgroundColor','none');


figure('Renderer', 'painters', 'Position', [50 50 400 300])
subplot('Position', [0.1 0.8 0.5 0.1])
scatter([1:length(n_cells)], ones(length(n_cells)), n_cells+1, 'k', 'filled')
ylim([0.9 1.1])
xlim([0.5 length(n_cells)+0.5])
axis off
title(['normalized CELL lineage proportions ' num2str(TargetAge) '-cell' ]);
set(gca,'fontsize', 7);
set(gca, 'FontName', 'Arial');
subplot('Position', [0.1 0.4 0.5 0.4])
h = heatmap(patlognorm2, 'MissingDataColor',[0.8 0.8 0.8]);
colormap(flipud(bone))
caxis([0 0.3])
h.XData = [pats];
h.YData = [1:NumClusters];
grid off
s = struct(h);
s.XAxis.TickLabelRotation = 45;
set(gca,'fontsize', 7);
%set(gca, 'FontName', 'Arial');
xlabel('cell ID')
ylabel('cluster ID')
%print([save_dir '/cellIDsHeatmapNorm ' num2str(TargetAge) ' cell'],'-dsvg');
exportgraphics(gcf, [save_dir '/cellIDsHeatmapNorm2 ' num2str(TargetAge) '-cell.pdf'],'ContentType','vector','BackgroundColor','none');


figure('Renderer', 'painters', 'Position', [50 50 400 300])
% subplot('Position', [0.1 0.8 0.5 0.1])
% scatter([1:length(n_cells)], ones(length(n_cells)), n_cells+1, 'k', 'filled')
% ylim([0.9 1.1])
% xlim([0.5 length(n_cells)+0.5])
% axis off
% set(gca,'fontsize', 7);
% set(gca, 'FontName', 'Arial');
subplot('Position', [0.1 0.4 0.5 0.4])

h = heatmap(patlogdbnorm, 'MissingDataColor',[0.8 0.8 0.8]);
colormap(flipud(bone))
caxis([0 0.3])
h.XData = [pats];
h.YData = [1:NumClusters];
grid off
s = struct(h);
s.XAxis.TickLabelRotation = 45;
set(gca,'fontsize', 7);
%set(gca, 'FontName', 'Arial');
title(['double normalized lineage cluster proportions ' num2str(TargetAge) '-cell' ])
xlabel('cell ID')
ylabel('cluster ID')
%print([save_dir '/cellIDsHeatmapNorm ' num2str(TargetAge) ' cell'],'-dsvg');
exportgraphics(gcf, [save_dir '/cellIDsHeatmapDBNorm ' num2str(TargetAge) '-cell.pdf'],'ContentType','vector','BackgroundColor','none');


% figure('Renderer', 'painters', 'Position', [300 300 400 300]);
% bh = bar(patlognorm,'stacked');
% colors = parula(numel(bh));
% for p = 1:numel(bh)
%     bh(p).FaceColor = colors(p,:);
% end
% legend(pats)
% legend('Location','eastoutside')
% set(gca, 'XTickLabel', 1:NumClusters)
% set(gca,'XTickLabelRotation',45)
% set(gca,'fontsize', 7);
% %set(gca, 'FontName', 'Arial');
% xlabel('cluster')
% ylabel('frequency')
% box off
% title(['cluster lineage proportions ' num2str(TargetAge) ' cell'])
% %print([dir '/cellIDsStackedBargraph ' num2str(TargetAge) ' cell'],'-dsvg');
% exportgraphics(gcf, [save_dir '/cellIDsStackedBargraph ' num2str(TargetAge) '-cell.pdf'],'ContentType','vector','BackgroundColor','none');


%% plot euclidean clustergram of cell normalized data https://www.mathworks.com/help/bioinfo/ug/working-with-the-clustergram-function.html

cg_s = clustergram(patlognorm2', 'RowLabels', pats,...
                               'Cluster', 'Column',...
                               'ImputeFun', @knnimpute)

cg_s.ColumnPDist = 'euclidean';
colormap(flipud(pink))
set(0,'ShowHiddenHandles','on')
% Get all handles from root
allhnds = get(0,'Children');
% Find hearmap axis and change the font size
h = findall(allhnds, 'Tag', 'HeatMapAxes');
set(h, 'FontSize', 7)
set(h, 'Clim', [0,0.3])
set(gcf, 'Position', [250 250 400 300])

exportgraphics(gcf, [save_dir '/cellIDsclustergram ' num2str(TargetAge) '-cell.pdf'],'ContentType','vector','BackgroundColor','none');


%% plot euclidean clustergram of cluster normalized https://www.mathworks.com/help/bioinfo/ug/working-with-the-clustergram-function.html

cg_s = clustergram(patlognorm', 'RowLabels', pats,...
                               'Cluster', 'Row',...
                               'ImputeFun', @knnimpute)

cg_s.ColumnPDist = 'euclidean';
colormap(flipud(pink))
set(0,'ShowHiddenHandles','on')
% Get all handles from root
allhnds = get(0,'Children');
% Find hearmap axis and change the font size
h = findall(allhnds, 'Tag', 'HeatMapAxes');
set(h, 'FontSize', 7)
set(h, 'Clim', [0,0.3])
set(gcf, 'Position', [250 250 400 300])

exportgraphics(gcf, [save_dir '/clus_cellIDsclustergram ' num2str(TargetAge) '-cell.pdf'],'ContentType','vector','BackgroundColor','none');
%% plot euclidean clustergram of double normalized https://www.mathworks.com/help/bioinfo/ug/working-with-the-clustergram-function.html
% 
% cg_s = clustergram(patlogdbnorm', 'RowLabels', pats,...
%                                'Cluster', 'all',...
%                                'ImputeFun', @knnimpute)
% 
% cg_s.ColumnPDist = 'euclidean';
% colormap(flipud(pink))
% set(0,'ShowHiddenHandles','on')
% % Get all handles from root
% allhnds = get(0,'Children');
% % Find hearmap axis and change the font size
% h = findall(allhnds, 'Tag', 'HeatMapAxes');
% set(h, 'FontSize', 7)
% set(h, 'Clim', [0,0.3])
% set(gcf, 'Position', [250 250 400 300])
% 
% exportgraphics(gcf, [save_dir '/db_cellIDsclustergram ' num2str(TargetAge) '-cell.pdf'],'ContentType','vector','BackgroundColor','none');

%% plot cosine clustergram https://www.mathworks.com/help/bioinfo/ug/working-with-the-clustergram-function.html

% cg_s = clustergram(patlognorm2', 'RowLabels', pats,...
%                                'Cluster', 'Column',...
%                                'ImputeFun', @knnimpute, 'ColumnPDist', 'cosine')
% 
% cg_s.ColumnPDist = 'cosine';
% cg_s.Colormap = redbluecmap;
% set(0,'ShowHiddenHandles','on')
% % Get all handles from root
% allhnds = get(0,'Children');
% % Find hearmap axis and change the font size
% h = findall(allhnds, 'Tag', 'HeatMapAxes');
% set(h, 'FontSize', 7)
% 
% set(gcf, 'Position', [250 250 400 300])
% 
% exportgraphics(gcf, [save_dir '/cellIDsclustergram_cosine ' num2str(TargetAge) '-cell.pdf'],'ContentType','vector','BackgroundColor','none');

%% kstest for distribution differences and pearson corr between cells
mkdir('stats')
lineage = save_dir;
save_dir = 'stats';

ks_pval = nan(length(pats), length(pats));


for i = 1:length(pats)
    for j = 1:length(pats)
        if i~=j
            [h,p] = kstest2(patlognorm2(:,i), patlognorm2(:,j));
            ks_pval(i,j) = p;
            if p<0.05
                disp([ string(pats(1,i)) ' cell ' 'and ' string(pats(1,j)) 'cell have different distributions of clusters p =' num2str(p)])
            end
        end
    end
end

figure('Renderer', 'painters', 'Position', [50 50 500 500])
h = heatmap(ks_pval, 'MissingDataColor',[0.8 0.8 0.8]);
colormap(flipud(bone))
caxis([0 0.1])
h.XData = [pats];
h.YData = [pats];
grid off
s = struct(h);
s.XAxis.TickLabelRotation = 45;
set(gca,'fontsize', 7);
title(['KS test pvalues ' num2str(TargetAge) '-cell' ])
xlabel('cell ID')
ylabel('cell ID')
%print([save_dir '/cellIDsHeatmapNorm ' num2str(TargetAge) ' cell'],'-dsvg');
exportgraphics(gcf, [save_dir '/' lineage '_cellIDsKSpval ' num2str(TargetAge) '-cell.pdf'],'ContentType','vector','BackgroundColor','none');

rho = corr(patlognorm2);
figure('Renderer', 'painters', 'Position', [50 50 500 500])
h = heatmap(rho, 'MissingDataColor',[0.8 0.8 0.8]);
colormap redbluecmap
%colormap(flipud(bone))
caxis([-1 1])
h.XData = [pats];
h.YData = [pats];
grid off
s = struct(h);
s.XAxis.TickLabelRotation = 45;
set(gca,'fontsize', 7);
title([lineage '_pearson corr ' num2str(TargetAge) '-cell' ])
xlabel('cell ID')
ylabel('cell ID')
exportgraphics(gcf, [save_dir '/' lineage '_cellIDsPearsonCorr ' num2str(TargetAge) '-cell.pdf'],'ContentType','vector','BackgroundColor','none');

%% compare sisters
sisters_corr = {};
nonsisters_corr = {};

for i = 1:length(pats)
    for j = 1:length(pats)
        if i~=j
            str1 = char(pats(1,i));
            str2 = char(pats(1,j));
            % 1. Check if they have the same length
            % 2. Check if all characters except the last are the same
            % 3. Check if the last characters are different
            isDiffByLast = (strlength(str1) == strlength(str2)) && ...
                strncmp(str1, str2, length(str1)-1) && ...
                (str1(end) ~= str2(end));
            if isDiffByLast ==1
                sisters_corr = [sisters_corr rho(i,j)]; %save pearson corr coef in sisters corr list
            else
                nonsisters_corr = [nonsisters_corr rho(i,j)]; %save pearson corr coef in nonsisters corr list
            end
        end
    end
end

% plot as boxplot - On each box, the central mark indicates the median,...
% and the bottom and top edges of the box indicate the 25th and 75th percentiles, respectively. ...
% The whiskers extend to the most extreme data points not considered outliers, and the outliers are plotted individually using the '+' marker symbol.
AllCorr{1,1} = cell2mat([sisters_corr]); 
AllCorr{1,2} = cell2mat([nonsisters_corr]);
% 1. Find the maximum length
maxlen = max(cellfun(@length, AllCorr));

% 2. Pad shorter vectors with NaN
C_padded = cellfun(@(V) [V, NaN(1, maxlen - length(V))], AllCorr, 'UniformOutput', false);

% 3. Optional: Convert to matrix
resultMatrix = cell2mat(C_padded'); 
meanCorr = nanmean(resultMatrix');
[h,p] = ttest2(AllCorr{1,1},AllCorr{1,2})


figure('Renderer', 'painters', 'Position', [300 300 150 150])
boxchart(resultMatrix','Notch','on');
hold on
plot(meanCorr,'-o')
set(gca,'fontsize', 7);
ylim([-1 1])
% xlabel('sister vs. non-sister cells')
% ylabel('Pearson correlation')
title([ lineage '-descendants p= ' num2str(p)])

exportgraphics(gcf, [save_dir '/' lineage '_cell_Pcorr_sisters.pdf'],'ContentType','vector','BackgroundColor','none');