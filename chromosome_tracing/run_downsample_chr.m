clear all
close all
save_path = 'downsampled_matrices/';
celltype = 'ABa'; 

load('AllChromosomes560_647_merge.mat');

TotalTADNum = 28;

%%
n = 0;
for i = 1:length(Chr)
    if startsWith(Chr(i).CellID, celltype)
        n = n+1;
        Chr_celltype(n)= Chr(i);
    end
end
clear n
length(Chr_celltype)

n = 115; %num to downsample randomly to
ChrDS = datasample(Chr_celltype, n);


STD = zeros(TotalTADNum,TotalTADNum);
Median = zeros(TotalTADNum,TotalTADNum);

for i = 1:TotalTADNum
    for j = 1:TotalTADNum
        DisList = [];
        if i ~= j
            for k = 1:length(ChrDS) %downsampled this value
                if ChrDS(k).r(i) == 1 && ChrDS(k).r(j) == 1
                   
                        DisList = [DisList ((ChrDS(k).x(i)-ChrDS(k).x(j))^2+(ChrDS(k).y(i)-ChrDS(k).y(j))^2+(ChrDS(k).z(i)-ChrDS(k).z(j))^2)^0.5];
                    
                end
            end
            
            STD(i,j) = std(DisList);
            Median(i,j) = median(DisList);
        end
    end
end


STD_total = zeros(TotalTADNum,TotalTADNum);
Median_total = zeros(TotalTADNum,TotalTADNum);

for i = 1:TotalTADNum
    for j = 1:TotalTADNum
        DisList2 = [];
        if i ~= j
            for k = 1:length(Chr_celltype)
                if Chr_celltype(k).r(i) == 1 && Chr_celltype(k).r(j) == 1
                   
                        DisList2 = [DisList2 ((Chr_celltype(k).x(i)-Chr_celltype(k).x(j))^2+(Chr_celltype(k).y(i)-Chr_celltype(k).y(j))^2+(Chr_celltype(k).z(i)-Chr_celltype(k).z(j))^2)^0.5];
                    
                end
            end
            
            STD_total(i,j) = std(DisList2);
            Median_total(i,j) = median(DisList2);
        end
    end
end

%%
figure('Renderer', 'painters', 'Position', [50 50 500 500])

ax1 = subplot(2,2,3);
hm1=heatmap(Median, 'MissingDataColor', [0.8 0.8 0.8]);
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
title([celltype '-descendants n = ' num2str(n) '/' num2str(length(Chr_celltype))])


ax2 = subplot(2,2,4);
hm2=heatmap(STD, 'MissingDataColor', [0.8 0.8 0.8]);
colormap(hm2,"parula")
grid off
caxis([0 1])
colorbar off
%mark the chromosome boundary
origState = warning('query', 'MATLAB:structOnObject');
cleanup = onCleanup(@()warning(origState));
warning('off','MATLAB:structOnObject')
S = struct(hm2); % Undocumented
axB = S.Axes;    % Undocumented
clear('cleanup')
S.XAxis.TickLabelRotation = 45;
set(gca,'fontsize', 7);
%set(gca, 'FontName', 'Arial');
mark = [TotalTADNum/2];
xline(axB, [mark+.5, mark+.5], 'k-');
yline(axB, [mark+.5, mark+.5], 'k-');
Labels = 1:TotalTADNum;
CustomLabels = string(Labels);
CustomLabels(mod(Labels,4) ~= 0) = " ";%show every 4 tick labels
hm2.XDisplayLabels = CustomLabels;
hm2.YDisplayLabels = CustomLabels;


%exportgraphics(gcf, [save_path 'DS_' num2str(n) '_median_std_' celltype '.pdf'],'ContentType','vector','BackgroundColor','none');

%figure('Renderer', 'painters', 'Position', [50 50 600 250])

ax3 = subplot(2,2,1);
hm1=heatmap(Median_total, 'MissingDataColor', [0.8 0.8 0.8]);
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
title([celltype '-descendants n = ' num2str(length(Chr_celltype))])


ax4 = subplot(2,2,2);
hm2=heatmap(STD_total, 'MissingDataColor', [0.8 0.8 0.8]);
colormap(hm2,"parula")
grid off
caxis([0 1])
colorbar off
%mark the chromosome boundary
origState = warning('query', 'MATLAB:structOnObject');
cleanup = onCleanup(@()warning(origState));
warning('off','MATLAB:structOnObject')
S = struct(hm2); % Undocumented
axB = S.Axes;    % Undocumented
clear('cleanup')
S.XAxis.TickLabelRotation = 45;
set(gca,'fontsize', 7);
%set(gca, 'FontName', 'Arial');
mark = [TotalTADNum/2];
xline(axB, [mark+.5, mark+.5], 'k-');
yline(axB, [mark+.5, mark+.5], 'k-');
Labels = 1:TotalTADNum;
CustomLabels = string(Labels);
CustomLabels(mod(Labels,4) ~= 0) = " ";%show every 4 tick labels
hm2.XDisplayLabels = CustomLabels;
hm2.YDisplayLabels = CustomLabels;


exportgraphics(gcf, [save_path 'total_ds_' num2str(n) '_median_std_' celltype '.pdf'],'ContentType','vector','BackgroundColor','none');

