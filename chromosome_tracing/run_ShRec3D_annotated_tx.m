%ShRec3d code to model average 3D structure from a contact map from Julien Mozziconacci (https://sites.google.com/site/julienmozziconacci/home/softwares?pli=1&authuser=0)

clear all
close all

TotalTADNum = 28;
save_path = 'sub0res06_sub1res08/';

%%%%specify block regions by cluster 11x28 matrix?
%dots = [8:24];

pat = 'D';
NumClusters = 11;

%load zygotic transcript signal
TX = load(['/scicore/home/mangos/finatu15/LineageLabeling/Analysis/Cole2024_GeneExpression_Reanalysis/TXsignal_lineages_zygotic/log10MeanTxsignal_' pat '.csv']);

for h=1:NumClusters
    
load([save_path 'ChrClus' num2str(h) '.mat'])
Median = zeros(TotalTADNum,TotalTADNum);

for i = 1:TotalTADNum
    for j = 1:TotalTADNum
        DisList = [];
        for k = 1:length(ChrClus)
            if ChrClus(k).r(i) == 1 && ChrClus(k).r(j) == 1
                if ((ChrClus(k).x(i)-ChrClus(k).x(j))^2+(ChrClus(k).y(i)-ChrClus(k).y(j))^2+(ChrClus(k).z(i)-ChrClus(k).z(j))^2)^0.5 <=6
                    DisList = [DisList ((ChrClus(k).x(i)-ChrClus(k).x(j))^2+(ChrClus(k).y(i)-ChrClus(k).y(j))^2+(ChrClus(k).z(i)-ChrClus(k).z(j))^2)^0.5];
                end
            end
        end
        Median(i,j) = median(DisList);
        
    end
end

Data=Median;

%% plot distance by tx signal product
Median_distances = [];
TX_product = [];
n = 0;
for i = 1:TotalTADNum
    for j = 1:TotalTADNum
        if i ~= j
            n = n+1;
            Median_distances(1,n) = Median(i,j);
            TX_product(1,n) = TX(i,1)*TX(j,1);
        end
    end
end
figure('Renderer', 'painters', 'Position', [300 300 400 250])
mdl = fitlm(Median_distances, TX_product); %fit linear regression model
stats = anova(mdl, 'summary');
p = plot(mdl);
p(1).Color = 'k'; p(1).Marker = 'o';
p(2).Color = 'k'; p(2).LineWidth = 1;
p(3).Color = 'k'; p(3).LineWidth = 1;
xlabel('pairwise distance')
ylabel('TX')
title(['cluster ' num2str(h)])
set(gca,'fontsize', 7);


%%

n = size(Data, 1);	
for k=1:n
i2k = repmat(Data(:,k), 1, n);
k2j = repmat(Data(k,:), n, 1);
Data = min(Data, i2k+k2j);
end   
n=length(Data);
center=zeros(n,1);
for i=1:n
    for j=1:n
    center(i)=center(i)+Data(i,j)^2-1/n*(Data(j,j:n)*Data(j:n,j));
    end
end
center(center<0)=0;
center=sqrt(center/n);
distmat=1/2*(repmat(center.*center,1,n)+repmat((center.*center)',n,1)-Data.*Data);
[V,D]=eigs(distmat,4);
lambda=diag(D);
XYZ=[V(:,1)*sqrt(lambda(1)) V(:,2)*sqrt(lambda(2)) V(:,3)*sqrt(lambda(3))];
scale=100/max(max(XYZ));
XYZ=XYZ*scale;

%% plot matrix and 3D model annotated with TX signal

TotalRegionNum = 14; %per chr
j = spring(TotalRegionNum);
j2 = winter(TotalRegionNum);
j3 = TX(1:14,:);%copper(TotalRegionNum);
j4 = TX(15:28,:);

xyzchrI = XYZ(1:14,:);
points = fnplt(cscvn(XYZ(1:14,:)'));

figure('Renderer', 'painters', 'Position', [10 10 350 650])
tiledlayout(7,1);
ax1 = nexttile;
%colormap(ax1,"redbluecmap")
heatmap(TX')
colorbar off
% plot3(points(1,:),points(2,:),points(3,:),'Color', 'm', 'LineWidth',2)
% hold on
% scatter3(xyzchrI(:,1),xyzchrI(:,2),xyzchrI(:,3), 100, j, 'filled')
% clear points
% xyzchrII = XYZ(15:28,:);
% points = fnplt(cscvn(XYZ(15:28,:)'));
% plot3(points(1,:),points(2,:),points(3,:),'Color', 'c', 'LineWidth',2)
% scatter3(xyzchrII(:,1),xyzchrII(:,2),xyzchrII(:,3), 100, j2, 'filled')
% hold off
%axis off
grid off
set(gca,'fontsize', 7);
title(['TX signal ' pat])

ax2 = nexttile([3 1]);
%colormap(ax2,"hot")
heatmap(Median)
clim([0 2])
grid off
set(gca,'fontsize', 7);

ax3 = nexttile([3 1]);
xyzchrI = XYZ(1:14,:);
points = fnplt(cscvn(XYZ(1:14,:)'));
plot3(points(1,:),points(2,:),points(3,:),'Color', 'm', 'LineWidth',2)
hold on
scatter3(xyzchrI(:,1),xyzchrI(:,2),xyzchrI(:,3), 200, j, 'LineWidth', 2)
scatter3(xyzchrI(:,1),xyzchrI(:,2),xyzchrI(:,3), 100, j3, 'filled')
clear points
xyzchrII = XYZ(15:28,:);
points = fnplt(cscvn(XYZ(15:28,:)'));
plot3(points(1,:),points(2,:),points(3,:),'Color', 'c', 'LineWidth',2)
scatter3(xyzchrII(:,1),xyzchrII(:,2),xyzchrII(:,3), 200, j2, 'LineWidth', 2)
scatter3(xyzchrII(:,1),xyzchrII(:,2),xyzchrII(:,3), 100, j4, 'filled')
colormap redbluecmap
colorbar
% for q=1:length(dots)
%     scatter3(XYZ(dots(1,q),1),XYZ(dots(1,q),2),XYZ(dots(1,q),3), 700,'MarkerEdgeColor', 'k')
% end
grid off
axis off
set(gca,'fontsize', 7);
%set(gca, 'FontName', 'Arial');
title(['cluster ' num2str(h)])
hold off

exportgraphics(gcf, [save_path '3D_cluster_' num2str(h) '_TX_' pat '.pdf'],'ContentType','vector','BackgroundColor','none');


end
