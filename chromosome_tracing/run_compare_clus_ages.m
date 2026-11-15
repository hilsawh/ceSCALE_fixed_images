clear all 
close all
save_path = 'sub0res06_sub1res08/';
TotalTADNum = 28;
NumClusters = 11;
%%
AllAges = {};

for i = 1:NumClusters
    load([save_path 'ChrClus' num2str(i) '.mat']);
    AllAges{1, i} = [ChrClus.Age];
end


%% plot age distributions as boxplot - On each box, the central mark indicates the median,...
% and the bottom and top edges of the box indicate the 25th and 75th percentiles, respectively. ...
% The whiskers extend to the most extreme data points not considered outliers, and the outliers are plotted individually using the '+' marker symbol.

% 1. Find the maximum length
maxlen = max(cellfun(@length, AllAges));

% 2. Pad shorter vectors with NaN
C_padded = cellfun(@(V) [V, NaN(1, maxlen - length(V))], AllAges, 'UniformOutput', false);

% 3. Optional: Convert to matrix
resultMatrix = cell2mat(C_padded'); 
meanAge = nanmean(resultMatrix');

figure('Renderer', 'painters', 'Position', [300 300 250 250])
boxchart(resultMatrix','Notch','on');
hold on
plot(meanAge,'-o')
set(gca,'fontsize', 7);
xlabel('cluster')
ylabel('age')
exportgraphics(gcf, [save_path 'cluster ages.pdf'],'ContentType','vector','BackgroundColor','none');

%%
[p,tbl,stats] = kruskalwallis(resultMatrix');
c = multcompare(stats);
tbl2 = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);