%Run this after run_scaling_PCAscore
close all
clear all

load('score_chr1.mat')

TotalNumTADs = 14;
k = 2; % the number of clusters to identify
r = 2; % the number of PCs to use
[~,V_temp,D_temp] = spectralcluster(score(:,1:r),5)

score1 = score(:,1);
score2 = score(:,2);
str = string(1:TotalNumTADs);

[idx1,V,D] = spectralcluster(score(:,1:r),k);
figure
j = jet(TotalNumTADs);
scatter(score1,score2,200,j, 'filled');
hold on
%textscatter(score1+0.01,score2+0.01,str);
gscatter(score1,score2,idx1, 'kkk', 'osx', 20);
axis square
Range = max([max(max(score(:,1)))-0.2 0.2-min(min(score(:,1)))]);
xlim([0-Range 0+Range])
Range = max([max(max(score(:,2)))-0.2 0.2-min(min(score(:,2)))]);
ylim([0-Range 0+Range])
grid on
xlabel('PC1')
ylabel('PC2')
title('spectralcluster compartment calls')
axis square
PlotProp
print('-painters','spectralcluster compartment calls','-dsvg')

idx2 = kmeans(score(:,1:r),k);
figure
scatter(score1,score2,200,j, 'filled');
hold on
%textscatter(score1+0.01,score2+0.01,str);
gscatter(score1,score2,idx2, 'kkk','osx', 20);
axis square
Range = max([max(max(score(:,1)))-0.2 0.2-min(min(score(:,1)))]);
xlim([0-Range 0+Range])
Range = max([max(max(score(:,2)))-0.2 0.2-min(min(score(:,2)))]);
ylim([0-Range 0+Range])
grid on
xlabel('PC1')
ylabel('PC2')
title('kmeans compartment calls')
axis square
PlotProp
print('-painters','kmeans compartment calls','-dsvg')

figure('Renderer', 'painters', 'Position', [800 800 400 400])
subplot(3,1,1)
bar(score1, 'w')
Range = 0.2;
ylim([0-Range 0+Range])
ylabel('PC1')
PlotProp
subplot(3,1,2)
bar(score2, 'w')
Range = 0.2;
ylim([0-Range 0+Range])
ylabel('PC2')
PlotProp
subplot(3,1,3)
imagesc(idx1')
colormap gray
ylabel('call')
xlabel('TADs')
PlotProp
print('-painters','compiledPC12spectral_2','-dsvg')

figure('Renderer', 'painters', 'Position', [800 800 400 400])
imagesc(score')
colormap gray
colorbar
title('score')
ylabel('PC')
xlabel('TADs')
PlotProp
print('-painters','score','-dsvg')

% T = tsne(score);
% figure
% scatter(T(:,1),T(:,2),200,j, 'filled')
% hold on
% textscatter(T(:,1),T(:,2), str)
% gscatter(T(:,1),T(:,2),idx2, 'kkk','osx', 20);
% grid on
% title('tsne tads')
% PlotProp


% figure('Renderer', 'painters', 'Position', [800 800 400 400])
% subplot(3,1,1)
% imagesc(score1')
% colormap gray
% %caxis([-0.5 0.5])
% colorbar
% ylabel('PC1')
% PlotProp
% subplot(3,1,2)
% imagesc(score2')
% colormap gray
% %caxis([-0.5 0.5])
% colorbar
% ylabel('PC2')
% PlotProp
% subplot(3,1,3)
% imagesc(idx1')
% colormap gray
% colorbar
% ylabel('call')
% xlabel('TADs')
% PlotProp
% print('-painters','compiledPC12spectral_1','-dsvg')

