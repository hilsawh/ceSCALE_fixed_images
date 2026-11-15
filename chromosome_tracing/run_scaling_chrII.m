clear all
close all

chromosome = 'chrII';
PCAonNormDisMatrix = 0;
NormalizeTowardsMean = 0;

load('../AllChromosomes_clusterannotated.mat');

TotalTADNum = 14; 
TotalTADNumBoth = 28;
chrIstart = 1;
chrIstop = 14;
chrIIstart = 15;
chrIIstop = 28;
AgeStart = 4;
AgeStop = 200;
Mean = zeros(TotalTADNumBoth,TotalTADNumBoth);
Std = zeros(TotalTADNumBoth,TotalTADNumBoth);
Median = zeros(TotalTADNumBoth,TotalTADNumBoth);

for i = chrIstart:chrIIstop
    for j = chrIstart:chrIIstop
        DisList = [];
        for k = 1:length(Chr) 
            if Chr(k).r(i) == 1 && Chr(k).r(j) == 1 && Chr(k).Age <= AgeStop && Chr(k).TracesInTerritory == 1    %pick only well-separated chromosomes            
                    DisList = [DisList ((Chr(k).x(i)-Chr(k).x(j))^2+(Chr(k).y(i)-Chr(k).y(j))^2+(Chr(k).z(i)-Chr(k).z(j))^2)^0.5];                
            end
        end
        Std(i,j) = std(DisList);
        SEM(i,j) = std(DisList)/(length(DisList))^0.5;
        NofData(i,j) = length(DisList);
        DisListAll{i,j} = DisList;
        Median(i,j) = median(DisList);
    end
end


Std1 = Std(chrIstart:chrIstop, chrIstart:chrIstop);
SEM1 = SEM(chrIstart:chrIstop, chrIstart:chrIstop);
NofData1 = NofData(chrIstart:chrIstop, chrIstart:chrIstop);
DisListAll1 = DisListAll(chrIstart:chrIstop, chrIstart:chrIstop);
Median1 = Median(chrIstart:chrIstop, chrIstart:chrIstop);


Std2 = Std(chrIIstart:chrIIstop, chrIIstart:chrIIstop);
SEM2 = SEM(chrIIstart:chrIIstop, chrIIstart:chrIIstop);
NofData2 = NofData(chrIIstart:chrIIstop, chrIIstart:chrIIstop);
DisListAll2 = DisListAll(chrIIstart:chrIIstop, chrIIstart:chrIIstop);
Median2 = Median(chrIIstart:chrIIstop, chrIIstart:chrIIstop);

%% chrII scaling

Median = Median2;
Std = Std2;
SEM = SEM2;
Std = Std2;
NofData = NofData2;
DisListAll = DisListAll2;


fid = fopen('chrIIprobeDomainStartsAndEnds.txt', 'r');
tline = fgetl(fid);
i = 0;
while ischar(tline)
    i = i+1;
    C = strsplit(tline);
    DomainStart(i) = str2num(C{2});
    DomainEnd(i) = str2num(C{3});
    DomainCenter(i) = (DomainStart(i)+DomainEnd(i))/2;
    tline = fgetl(fid);
end
fclose(fid);

DomainCenter = DomainCenter./1000000; %convert from bp to Mb

GenomicDis = zeros(TotalTADNum,TotalTADNum);
GenomicDis_pxl = [];
Dis_pxl = [];
Std_pxl = [];
CV_pxl = [];

for i = 1:TotalTADNum
    for j = 1:TotalTADNum
        GenomicDis(i,j) = abs(DomainCenter(i)-DomainCenter(j));
        if i~=j
            GenomicDis_pxl = [GenomicDis_pxl GenomicDis(i,j)];
            Dis_pxl = [Dis_pxl Median(i,j)];
            Std_pxl = [Std_pxl Std(i,j)];            
        end
    end
end

%% scaling plots

f = figure('Renderer', 'painters', 'Position', [10 10 250 250]);
f.Name = '1';
set(gca,'fontsize', 8);
heatmap(GenomicDis)
colormap hot
colorbar
title('genomic distance - Mb');

f = figure('Renderer', 'painters', 'Position', [10 10 250 250]);
f.Name = '2';
plot(GenomicDis_pxl, Dis_pxl, 'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w');
xlabel('genomic distance (Mb)');
ylabel('mean spatial distance (um)');
R = corrcoef(GenomicDis_pxl, Dis_pxl);

% fit to the fractal globule model
x = log(GenomicDis_pxl)';
y = log(Dis_pxl)'-x/3;
X = [ones(size(x))];
b = regress(y,X);
x = GenomicDis_pxl';
y_fit = exp(b)*x.^(1/3);

f = figure('Renderer', 'painters', 'Position', [10 10 250 250]);
f.Name = '3';
loglog(GenomicDis_pxl, Dis_pxl, 'o', 'MarkerSize', 5);
xlabel('log(genomic distance)');
ylabel('log(median spatial distance)');hold on
[x, IX] = sort(x);
y_fit = y_fit(IX);
plot(x, y_fit, 'k--', 'LineWidth', 2);
hold off

x = log(GenomicDis_pxl)';
y = log(Dis_pxl)';
X = [ones(size(x)), x];
[b, bint] = regress(y,X);
%save(['scaling.mat'], 'b');
title(['Slope = ' num2str(b(2)) '+-' num2str(bint(2,2)-b(2))])
x = GenomicDis_pxl';
y_fit = exp(b(1))*x.^b(2);

figure(3)
title(['Slope = ' num2str(b(2)) '+-' num2str(bint(2,2)-b(2))])
hold on
[xp, IX] = sort(x);
y_fitp = y_fit(IX);
plot(xp, y_fitp, 'k-', 'LineWidth', 2);
hold off
xlim([0 15])
ylim([0.4 2.5])


%% Normalize median distance to genomic expected
[xData, yData] = prepareCurveData( GenomicDis_pxl, Dis_pxl );

% Set up fittype and options.
ft = fittype( 'power1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.939816085593886 0.24545186726114];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
f = figure('Renderer', 'painters', 'Position', [10 10 250 250]);
f.Name = '8';
set(gca,'fontsize', 8);
h = plot(fitresult, xData, yData, 'predobs');
% legend( h, 'Dis vs. GenomicDis', 'distance scaling', 'Lower bounds (distance scaling)', 'Upper bounds (distance scaling)', 'Location', 'SouthEast' );
% Label axes
xlabel('genomic distance')
ylabel('spatial distance')
grid off
xlim([0 15])
axis square
exportgraphics(gcf, ['Scaling_' chromosome '.pdf'],'ContentType','vector','BackgroundColor','none');


fitresult
fitcoef=coeffvalues(fitresult);
title(['y = ' num2str(fitcoef(1,1)) 'x^' num2str(fitcoef(1,2))]);


x = GenomicDis_pxl';
y_fit2 = fitcoef(1,1)*x.^fitcoef(1,2);

Median_adjust = zeros(TotalTADNum,TotalTADNum);
n = 0;
if NormalizeTowardsMean == 0
    for i = 1:TotalTADNum
        for j = 1:TotalTADNum
            if i~=j
                n = n+1;
                Median_adjust(i,j) = Median(i,j)/y_fit2(n);
            else
                Median_adjust(i,j) = 1;
            end
        end
    end
else
    %calculate the running average of spatial distance, 
    %use the 40 closest genomic distance entries for each data point.
    for i = 1:TotalTADNum
        for j = 1:TotalTADNum
            if i~=j
                n = n+1;
                GenomicDis_diff = abs(GenomicDis_pxl - GenomicDis_pxl(n));
                [GenomicDis_sort, IX] = sort(GenomicDis_diff);
                Dis_sort = Dis_pxl(IX);
                Dis_sort_mean(n) = mean(Dis_sort(1:40));
                Median_adjust(i,j) = Median(i,j)/Dis_sort_mean(n);
            else
                Median_adjust(i,j) = 1;
            end
        end
    end
    figure
    plot(GenomicDis_pxl, Dis_pxl, '.', 'MarkerSize', 20);
    xlabel('pairwise genomic distance (Mb)');
    ylabel('pairwise spatial distance (um)');

    [GenomicDis_sort, IX] = sort(GenomicDis_pxl);
    hold on
    plot(GenomicDis_sort, Dis_sort_mean(IX), 'k-', 'LineWidth', 5);
    hold off
    xlim([0 21])
    axis square
end

f = figure('Renderer', 'painters', 'Position', [10 10 250 250]);
f.Name = '5';
set(gca,'fontsize', 8);
heatmap(Median_adjust)
Range = max([max(max(Median_adjust))-1 1-min(min(Median_adjust))]);
clim([0.8 1.2])
colormap redbluecmap
colorbar
grid off
title('nomalized spatial distance');
exportgraphics(gcf, ['normalized_distance' chromosome '.pdf'],'ContentType','vector','BackgroundColor','none');


% calculate correlation matrix
Pearson = zeros(TotalTADNum,TotalTADNum);

for i = 1:TotalTADNum
    for j = 1:TotalTADNum
        X = [];
        Y = [];
        for k = 1:TotalTADNum
            if k~=i && k~=j
                X = [X; Median_adjust(k,i)];
                Y = [Y; Median_adjust(k,j)];
            end
        end
        Pearson(i,j) = corr(X,Y, 'type', 'Pearson');
        
    end
end

f = figure('Renderer', 'painters', 'Position', [10 10 250 250]);
f.Name = '6';
set(gca,'fontsize', 8);
heatmap(Pearson)
colormap redbluecmap
colorbar
grid off
% Range = max([max(max(Pearson))-1 1-min(min(Pearson))]);
% clim([0-Range 0+Range])
title('Pearson correlation')
exportgraphics(gcf, ['Pearson_' chromosome '.pdf'],'ContentType','vector','BackgroundColor','none');


if PCAonNormDisMatrix == 0
    [coeff,score,latent,tsquared,explained,mu] = pca(Pearson);
else
	[coeff,score,latent,tsquared,explained,mu] = pca(Median_adjust);
end

explained;

f = figure('Renderer', 'painters', 'Position', [10 10 250 250]);
f.Name = '7';
set(gca,'fontsize', 8);
j = winter(TotalTADNum);
scatter3(score(:,1),score(:,2),score(:,3),100,j, 'filled');
labels = {'15' '16' '17' '18' '19' '20' '21' '22' '23' '24' '25' '26' '27' '28'};
text(score(:,1),score(:,2),score(:,3), labels, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
xlabel('PC 1')
ylabel('PC 2')
zlabel('PC 3')
Range = max([max(max(score(:,1)))-0.2 0.2-min(min(score(:,1)))]);
xlim([0-Range 0+Range])
Range = max([max(max(score(:,2)))-0.2 0.2-min(min(score(:,2)))]);
ylim([0-Range 0+Range])
colormap winter
colorbar
c = colorbar;
c.TickDirection = 'out';
c.Location = 'northoutside';
clim([1 TotalTADNum])
c.Ticks = [1:TotalTADNum];
c.Label.String = 'regions';
view(2)
exportgraphics(gcf, ['PCA_3D_' chromosome '.pdf'],'ContentType','vector','BackgroundColor','none');


f = figure('Renderer', 'painters', 'Position', [10 10 250 250]);
f.Name = '8';
set(gca,'fontsize', 8);
bar(explained)
axis square
title('PCA explained')
xlabel('PCs')
xlim([0 15])
ylabel('% variation explained')

exportgraphics(gcf, ['PCA_explained_' chromosome '.pdf'],'ContentType','vector','BackgroundColor','none');

%plot pc1,2,3 to compare to hi-c A/B

f = figure('Renderer', 'painters', 'Position', [10 10 250 250]);
f.Name = '9';
set(gca,'fontsize', 8);
subplot(3,1,1)
bar(score(:,1), 'w')
xlim([0 15])
xticklabels({});
ylabel('PC1')
box off
subplot(3,1,2)
bar(score(:,2), 'w')
xlim([0 15])
xticklabels({});
ylabel('PC2')
box off
subplot(3,1,3)
bar(score(:,3), 'w')
xlim([0 15])
ylabel('PC3')
box off

exportgraphics(gcf, ['PCA_compartments_' chromosome '.pdf'],'ContentType','vector','BackgroundColor','none');
