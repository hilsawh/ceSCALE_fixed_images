clear all
close all

load('Clustering/AllChromosomes_clusterannotated.mat');

TotalTADNum = 28;
NumClusters = 11;

agebins(:,1) = [0:25:75];
agebins(:,2) = [agebins(2,1)-1:25:100];

agebinnames = string(zeros(length(agebins), 1));
for i = 1:length(agebins)
    agebinnames(i,1) = string([num2str(agebins(i,1)) '-' num2str(agebins(i,2))]);
end

%% calculate median distance matrix over time
AllMedians = [];
All_n = [];

for p = 1:length(agebins)
    agestart = agebins(p,1);
    agestop = agebins(p,2);
    n = 0;
    for q = 1:length(Chr)
        if Chr(q).Age >= agestart && Chr(q).Age <= agestop
            n = n+1;
        end
    end
    All_n = [All_n n];
end
disp([num2str(sum(All_n)) ' total traces'])

for p = 1:length(agebins)
    Median = zeros(TotalTADNum, TotalTADNum);
    agestart = agebins(p,1);
    agestop = agebins(p,2);
    for i = 1:TotalTADNum
        for j = 1:TotalTADNum
            DisList = [];
            if i ~= j
                for q = 1:length(Chr)
                    if Chr(q).Age >= agestart && Chr(q).Age <= agestop
                        if Chr(q).r(i) == 1 && Chr(q).r(j) == 1
                            DisList = [DisList ((Chr(q).x(i)-Chr(q).x(j))^2+(Chr(q).y(i)-Chr(q).y(j))^2+(Chr(q).z(i)-Chr(q).z(j))^2)^0.5];
                        end
                    end
                end
                Median(i,j) = median(DisList);
            end
        end
    end
    AllMedians{1,p} = Median;
end


% plot matrices

figure('Renderer', 'painters', 'Position', [50 50 800 200])
for i = 1:length(agebins)
    ax = subplot(1, length(agebins), i);
    hm1=heatmap(AllMedians{1,i}, 'MissingDataColor', [0.8 0.8 0.8]);
    colormap(hm1, "hot")
    grid off
    caxis([0 2])
    colorbar off
    %mark the chromosome boundary
    origState = warning('query', 'MATLAB:structOnObject');
    cleanup = onCleanup(@()warning(origState));
    warning('off','MATLAB:structOnObject')
    S = struct(hm1); % Undocumented
    axA = S.Axes;    % Undocumented
    clear('cleanup')
    S.XAxis.TickLabelRotation = 45;
    set(gca,'fontsize', 7);
    % %set(gca, 'FontName', 'Arial');
    mark = [TotalTADNum/2];
    xline(axA, [mark+.5, mark+.5], 'k-');
    yline(axA, [mark+.5, mark+.5], 'k-');
    Labels = 1:TotalTADNum;
    CustomLabels = string(Labels);
    CustomLabels(mod(Labels,4) ~= 0) = " ";%show every 4 tick labels
    hm1.XDisplayLabels = CustomLabels;
    hm1.YDisplayLabels = CustomLabels;
    title([agebinnames(i,1) num2str(All_n(1,i))])
end
colorbar

exportgraphics(gcf, 'median_dist_over_time.pdf','ContentType','vector','BackgroundColor','none');

%% downsample traces in each age bin
AllMedians = [];
for p = 1:length(agebins)
    Median = zeros(TotalTADNum, TotalTADNum);
    agestart = agebins(p,1);
    agestop = agebins(p,2);
    
    n = 0;
    for q = 1:length(Chr)
        if Chr(q).Age >= agestart && Chr(q).Age <= agestop
            n = n+1;
            Chr_bin(n)= Chr(q);
        end
    end
    clear n
    
    n = 414; %num to downsample randomly to
    ChrDS = datasample(Chr_bin, n);
    for i = 1:TotalTADNum
        for j = 1:TotalTADNum
            DisList = [];
            if i ~= j
                for q = 1:length(ChrDS)
                    %if ChrDS(q).Age >= agestart && ChrDS(q).Age <= agestop
                        if ChrDS(q).r(i) == 1 && ChrDS(q).r(j) == 1
                            DisList = [DisList ((ChrDS(q).x(i)-ChrDS(q).x(j))^2+(ChrDS(q).y(i)-ChrDS(q).y(j))^2+(ChrDS(q).z(i)-ChrDS(q).z(j))^2)^0.5];
                        end
                    %end
                end
                Median(i,j) = median(DisList);
            end
        end
    end
    AllMedians{1,p} = Median;
    clear Chr_bin
end
% plot matrices

figure('Renderer', 'painters', 'Position', [50 50 800 200])
for i = 1:length(agebins)
    ax = subplot(1, length(agebins), i);
    hm1=heatmap(AllMedians{1,i}, 'MissingDataColor', [0.8 0.8 0.8]);
    colormap(hm1, "hot")
    grid off
    caxis([0 2])
    colorbar off
    %mark the chromosome boundary
    origState = warning('query', 'MATLAB:structOnObject');
    cleanup = onCleanup(@()warning(origState));
    warning('off','MATLAB:structOnObject')
    S = struct(hm1); % Undocumented
    axA = S.Axes;    % Undocumented
    clear('cleanup')
    S.XAxis.TickLabelRotation = 45;
    set(gca,'fontsize', 7);
    % %set(gca, 'FontName', 'Arial');
    mark = [TotalTADNum/2];
    xline(axA, [mark+.5, mark+.5], 'k-');
    yline(axA, [mark+.5, mark+.5], 'k-');
    Labels = 1:TotalTADNum;
    CustomLabels = string(Labels);
    CustomLabels(mod(Labels,4) ~= 0) = " ";%show every 4 tick labels
    hm1.XDisplayLabels = CustomLabels;
    hm1.YDisplayLabels = CustomLabels;
    title([agebinnames(i,1) ' cell'])
end
colorbar

exportgraphics(gcf, 'median_dist_over_time_downsampled.pdf','ContentType','vector','BackgroundColor','none');

% difference map oldest to youngest

Diff = AllMedians{1,length(agebins)}- AllMedians{1,1};

figure('Renderer', 'painters', 'Position', [50 50 200 200])
hm1=heatmap(Diff, 'MissingDataColor', [0.8 0.8 0.8]);
colormap(flipud(redbluecmap))
grid off
caxis([-1 1])
%colorbar off
%mark the chromosome boundary
origState = warning('query', 'MATLAB:structOnObject');
cleanup = onCleanup(@()warning(origState));
warning('off','MATLAB:structOnObject')
S = struct(hm1); % Undocumented
axA = S.Axes;    % Undocumented
clear('cleanup')
S.XAxis.TickLabelRotation = 45;
set(gca,'fontsize', 7);
% %set(gca, 'FontName', 'Arial');
mark = [TotalTADNum/2];
xline(axA, [mark+.5, mark+.5], 'k-');
yline(axA, [mark+.5, mark+.5], 'k-');
Labels = 1:TotalTADNum;
CustomLabels = string(Labels);
CustomLabels(mod(Labels,4) ~= 0) = " ";%show every 4 tick labels
hm1.XDisplayLabels = CustomLabels;
hm1.YDisplayLabels = CustomLabels;
title([agebinnames(1,1) agebinnames(length(agebins),1)])


exportgraphics(gcf, 'median_dist_downsampled_differencemap.pdf','ContentType','vector','BackgroundColor','none');
