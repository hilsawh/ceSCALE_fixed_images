clear all
close all
HiC = load('DCC_WT_HiC.txt');
NumBins = size(HiC,1);
BinCenters = 50000*[1:NumBins] - 25000; %for 50kb bins, change this for different bin sizes

TotalTADNum = 43;
load('DisListAll.mat');
load('MeanPairDistanceMatrix561.mat');


figure
heatmap(HiC)
grid off
caxis([0 4])
title('HiC map chrV full matrix')
colormap(flipud(hot))

%subset full HiC matrix to get contacts for probed regions

fid = fopen('StartsAndEndschrX_probe_design.txt', 'r');
tline = fgetl(fid);
i = 0;
while ischar(tline)
    i = i+1;
    C = strsplit(tline);
    DomainStart(i) = str2num(C{2});
    DomainEnd(i) = str2num(C{3});
    BinCenter(i) = (DomainStart(i)+DomainEnd(i))/2;
    tline = fgetl(fid);
end
fclose(fid);
for i = 1:length(DomainStart)
    BinIDs = [];
    for j = 1:NumBins
        if BinCenters(j)>DomainStart(i) && BinCenters(j)<DomainEnd(i)
            BinIDs = [BinIDs j];
        end
    end
    ChozenBinID(i,:) = [min(BinIDs) max(BinIDs)];
end

for i = 1:TotalTADNum
    for j = 1:TotalTADNum
        HiC_chozen(i, j) = sum(sum(HiC(ChozenBinID(i,1):ChozenBinID(i,2),ChozenBinID(j,1):ChozenBinID(j,2))));
        HiC_chozen(i, j) = HiC_chozen(i, j)/(ChozenBinID(i,2)-ChozenBinID(i,1)+1)/(ChozenBinID(j,2)-ChozenBinID(j,1)+1);
    end
end

figure
for i = 1:length(HiC_chozen)
    HiC_chozen(i, i) = 0;
end
for i = 1:length(HiC_chozen)
    HiC_chozen(i, i) = max(max(HiC_chozen));
end
imagesc(HiC_chozen)
colorbar
title('HiC map probed regions'); 
axis square
colormap(flipud(hot))

%% compare spatial distance matrix to HiC matrix for FISH probed regions

figure
imagesc(Mean)
colorbar
title('Mean pairwise distance');
axis square
colormap hot
caxis([0 2])
PlotProp

GenomicDis = zeros(TotalTADNum,TotalTADNum);
GenomicDis_pxl = [];
HiC_pxl = [];
Dis_pxl = [];

CloseProb_pxl = [];
CloseThresh = 0.5;
for i = 1:TotalTADNum
    for j = 1:TotalTADNum 
        GenomicDis(i,j) = abs(BinCenter(i)-BinCenter(j));
        if i~= j  
            GenomicDis_pxl = [GenomicDis_pxl GenomicDis(i,j)];
            HiC_pxl = [HiC_pxl HiC_chozen(i,j)];
            Dis_pxl = [Dis_pxl Mean(i,j)];
            CloseProb_pxl = [CloseProb_pxl length(find(DisListAll{i,j}<CloseThresh))/length(DisListAll{i,j})];
        end
    end
end

figure(30)
loglog(Dis_pxl, 1./HiC_pxl, '.','MarkerSize', 20);
xlabel('Mean pariwise distance (um)');
ylabel('Inversed HiC signal');
PlotProp
hold on
axis square
Ind = find(HiC_pxl>0);
x = log(Dis_pxl(Ind))';
y = log(1./HiC_pxl(Ind))';
X = [ones(size(x)), x];
[b, bint] = regress(y,X);
display(['Slope = ' num2str(b(2)) '+-' num2str(bint(2,2)-b(2))])
R = corrcoef(x, y);
display(['Correlation coefficient = ' num2str(R(2,1))])
x = Dis_pxl';
y_fit = exp(b(1))*x.^b(2);
loglog(x, y_fit, '-r','LineWidth', 5);
hold off
bint;

%% pca analysis on HiC for FISH probed regions only

Mean = 1./HiC_chozen;

fid = fopen('chrVprobeDomainStartsAndEnds.txt', 'r');
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


for i = 1:TotalTADNum
    for j = 1:TotalTADNum
        GenomicDis(i,j) = abs(DomainCenter(i)-DomainCenter(j));
        if i~=j
            GenomicDis_pxl = [GenomicDis_pxl GenomicDis(i,j)];
            Dis_pxl = [Dis_pxl Mean(i,j)];
        end
    end
end


figure
imagesc(GenomicDis)
colormap hot
colorbar
title('Genomic distance - Mb');
PlotProp
axis square

%scale the data

figure(300)
plot(GenomicDis_pxl, Dis_pxl, '.', 'MarkerSize', 20);
xlabel('Genomic distance');
ylabel('inverse CF');
R = corrcoef(GenomicDis_pxl, Dis_pxl);
% title(['Correlation coefficient = ' num2str(R(2,1))])

figure(400)
loglog(GenomicDis_pxl, Dis_pxl, '.', 'MarkerSize', 20);
xlabel('log(Genomic distance)');
ylabel('log(inverse CF)');
PlotProp
hold on


% fit to get new scaling factor
figure(400)
hold on
x = log(GenomicDis_pxl)';
y = log(Dis_pxl)';
X = [ones(size(x)), x];
[b, bint] = regress(y,X);

title(['Slope = ' num2str(b(2)) '+-' num2str(bint(2,2)-b(2))])

x = GenomicDis_pxl';
y_fit = exp(b(1))*x.^b(2);

figure(300)
title(['Slope = ' num2str(b(2)) '+-' num2str(bint(2,2)-b(2))])
hold on
[xp, IX] = sort(x);
y_fitp = y_fit(IX);
plot(xp, y_fitp, 'r-', 'LineWidth', 5);
hold off
xlim([0 21])
% ylim([0.4 1.8])
axis square
PlotProp
%print('InvCFScaling','-dsvg')

figure(400)
loglog(xp, y_fitp, 'r-', 'LineWidth', 5);
hold off
xlim([0 21])
axis square
PlotProp
%print('LogInvCFScaling','-dsvg')


Mean_adjust = Mean;
n = 0;

for i = 1:TotalTADNum
    for j = 1:TotalTADNum
        if i~=j
            n = n+1;
            Mean_adjust(i,j) = Mean(i,j)/y_fit(n);
        else
            Mean_adjust(i,j) = 1;
        end
    end
end

figure
heatmap(Mean_adjust)
RedBlue
colorbar
caxis([1-1 1+1])
title('Scaled HiC probed regions only')

[coeff,score,latent,tsquared,explained,mu] = pca(Mean_adjust);

score1 = score(:,1);
figure
bar(score1);
ylabel('PC -1');
xlabel('TADs');
hold on
A = find(score1<0);
bar(A, score1(A), 'r');
hold off
PlotProp
    
%% pca analysis for HiC full matrix

Mean = 1./HiC;

GenomicDis = zeros(NumBins,NumBins);
GenomicDis_pxl = [];
Dis_pxl = [];

HiCbinlocations= zeros(NumBins,2);
for i = 2:NumBins
    HiCbinlocations(i,1) = HiCbinlocations(i-1) + 50001;
end
for i= 1:NumBins
    HiCbinlocations(i,2) = HiCbinlocations(i) + 50000;
    BinCenter(i) = (HiCbinlocations(i,1)+HiCbinlocations(i,2))/2;
end


for i = 1:NumBins
    for j = 1:NumBins
        if i~=j
        GenomicDis(i,j) = abs(BinCenter(i)-BinCenter(j));
        end
        if i~=j
            GenomicDis_pxl = [GenomicDis_pxl GenomicDis(i,j)];
            Dis_pxl = [Dis_pxl Mean(i,j)];
        end
    end
end


figure
imagesc(GenomicDis)
colormap hot
colorbar
title('Genomic distance');
PlotProp
axis square

%scale the data

figure(30)
plot(GenomicDis_pxl, Dis_pxl, '.', 'MarkerSize', 20);
xlabel('Genomic distance');
ylabel('inverse CF');
R = corrcoef(GenomicDis_pxl, Dis_pxl);
% title(['Correlation coefficient = ' num2str(R(2,1))])

figure(40)
loglog(GenomicDis_pxl, Dis_pxl, '.', 'MarkerSize', 20);
xlabel('log(Genomic distance)');
ylabel('log(inverse CF)');
PlotProp
hold on

% fit to get scaling factor
[xData, yData] = prepareCurveData( GenomicDis_pxl, Dis_pxl );

% Set up fittype and options.
ft = fittype( 'power1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.939816085593886 0.24545186726114];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'cf scaling' );
h = plot( fitresult, xData, yData, 'predobs' );
legend( h, 'invCF vs. GenomicDis', 'invCF scaling', 'Lower bounds (invCF scaling)', 'Upper bounds (invCF scaling)', 'Location', 'SouthEast' );
% Label axes
xlabel GenomicDis
ylabel invCF
grid off
axis square
PlotProp
fitresult
fitcoef=coeffvalues(fitresult);
title(['y = ' num2str(fitcoef(1,1)) 'x^' num2str(fitcoef(1,2))]);

Mean_adjust = zeros(NumBins,NumBins);
n=0;
for i = 1:NumBins
    for j = 1:NumBins
        if i~=j
            n=n+1;
            Mean_adjust(i,j) = Mean(i,j)/fitresult(GenomicDis(i,j));
        else
            Mean_adjust(i,j) = 1;
        end
    end
end

Mean_adjust(isinf(Mean_adjust)) = nan;

figure
heatmap(Mean_adjust)
grid off
RedBlue
caxis([1-2 1+2])
colorbar
title('Nomalized inverse CF');

Mean_adjust(isinf(Mean_adjust)| isnan(Mean_adjust)) = 1;

[coeff,score,latent,tsquared,explained,mu] = pca(Mean_adjust);
score1 = score(:,1);

figure
bar(score1);
ylabel('PC -1');
xlabel('bins');
hold on
A = find(score1<0);
bar(A, score1(A), 'r');
hold off
%ylim([-1 1])
xlim([0 NumBins+1])
xticks([1:10:NumBins])
xtickangle(45)
PlotProp
%print('PC1','-dsvg')






