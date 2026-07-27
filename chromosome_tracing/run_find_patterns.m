%% find patterns of chromosome conformation in different lineages

clear all
%close all

AgeStart = 2; %change these to specify embryo age
AgeStop = 150;
AgeRange = ([ num2str(AgeStart) 'to' num2str(AgeStop) 'cell']);
save_path = ['chr_matrix_' AgeRange];
load([save_path '/AllChromosomes560_647_merge.mat'])


pat = ['P2'];

newChrName = ['chr_' pat];
n = 0;
n_other = 0;
for i = 1:length(Chr)
    if strcmp(Chr(i).CellID, pat) %for exact pat
    %if startsWith(Chr(i).CellID, pat) %for all in lineage after pat
        n = n+1;
        struct(n)=Chr(i);
    else
        n_other = n_other+1;
        struct_other(n_other)=Chr(i);
    end
end
disp(['number of ' pat ' cells is ' num2str(n)])
disp(['number of other cells is ' num2str(n_other)])

TotalTADNumBoth = 28;
chrIstart = 1;
chrIstop = 14;
chrIIstart = 15;
chrIIstop = 28;

%% plot pat matrix

Mean = zeros(TotalTADNumBoth,TotalTADNumBoth);
Std = zeros(TotalTADNumBoth,TotalTADNumBoth);
Median = zeros(TotalTADNumBoth,TotalTADNumBoth);

for i = chrIstart:chrIIstop
    for j = chrIstart:chrIIstop
        if i ~= j
        DisList = [];
        
        for k = 1:length(struct) 
            if struct(k).r(i) == 1 && struct(k).r(j) == 1 %&& struct(k).TracesInTerritory == 1    %pick only well-separated chromosomes            
                    DisList = [DisList ((struct(k).x(i)-struct(k).x(j))^2+(struct(k).y(i)-struct(k).y(j))^2+(struct(k).z(i)-struct(k).z(j))^2)^0.5];                
            end
        end
        Mean(i,j) = mean(DisList);
        Std(i,j) = std(DisList);
        SEM(i,j) = std(DisList)/(length(DisList))^0.5;
        NofData(i,j) = length(DisList);
        DisListAll{i,j} = DisList;
        Median(i,j) = median(DisList);
        end
    end
end
save([save_path '/AllChromosomes560_647_merge_' pat '.mat'],'struct')

clear struct

figure('Renderer', 'painters', 'Position', [600 600 275 250])
hm=heatmap(Median, 'MissingDataColor', [0.8 0.8 0.8]);
colormap hot
grid off
caxis([0 2])
colorbar off
%mark the chromosome boundary
origState = warning('query', 'MATLAB:structOnObject');
cleanup = onCleanup(@()warning(origState));
warning('off','MATLAB:structOnObject')
S = struct(hm); % Undocumented
ax = S.Axes;    % Undocumented
clear('cleanup')
S.XAxis.TickLabelRotation = 45;
set(gca,'fontsize', 7);
%set(gca, 'FontName', 'Arial');
mark = [chrIstop];
xline(ax, [mark+.5, mark+.5], 'k-');
yline(ax, [mark+.5, mark+.5], 'k-');
title(['median spatial distance ' pat ' cells n=' num2str(n)]);
exportgraphics(gcf, [save_path '/median spatial distance' AgeRange '_' newChrName '.pdf'],'ContentType','vector','BackgroundColor','none');

figure('Renderer', 'painters', 'Position', [600 600 275 250])
hm=heatmap(Std, 'MissingDataColor', [0.8 0.8 0.8]);
colormap parula %abyss
grid off
caxis([0 1])
colorbar off
%mark the chromosome boundary
origState = warning('query', 'MATLAB:structOnObject');
cleanup = onCleanup(@()warning(origState));
warning('off','MATLAB:structOnObject')
S = struct(hm); % Undocumented
ax = S.Axes;    % Undocumented
clear('cleanup')
S.XAxis.TickLabelRotation = 45;
set(gca,'fontsize', 7);
%set(gca, 'FontName', 'Arial');
mark = [chrIstop];
xline(ax, [mark+.5, mark+.5], 'k-');
yline(ax, [mark+.5, mark+.5], 'k-');
title(['std ' pat ' cells'])
exportgraphics(gcf, [save_path '/std spatial distance' AgeRange '_' newChrName '.pdf'],'ContentType','vector','BackgroundColor','none');
% 
% figure
% heatmap(NofData)
% colormap gray
% title('number of measurements');
% exportgraphics(gcf, [save_path '/num measurements' AgeRange '_' newChrName '.pdf'],'ContentType','vector','BackgroundColor','none');

%% plot other cells matrix
% 
% Mean = zeros(TotalTADNumBoth,TotalTADNumBoth);
% Std = zeros(TotalTADNumBoth,TotalTADNumBoth);
% Median = zeros(TotalTADNumBoth,TotalTADNumBoth);
% 
% for i = chrIstart:chrIIstop
%     for j = chrIstart:chrIIstop
%         if i ~= j
%         DisList = [];
% 
%         for k = 1:length(struct_other) 
%             if struct_other(k).r(i) == 1 && struct_other(k).r(j) == 1 %&& struct_other(k).TracesInTerritory == 1    %pick only well-separated chromosomes            
%                     DisList = [DisList ((struct_other(k).x(i)-struct_other(k).x(j))^2+(struct_other(k).y(i)-struct_other(k).y(j))^2+(struct_other(k).z(i)-struct_other(k).z(j))^2)^0.5];                
%             end
%         end
%         Mean(i,j) = mean(DisList);
%         Std(i,j) = std(DisList);
%         SEM(i,j) = std(DisList)/(length(DisList))^0.5;
%         NofData(i,j) = length(DisList);
%         DisListAll{i,j} = DisList;
%         Median(i,j) = median(DisList);
%         end
%     end
% end
% % save([save_path '/AllChromosomes560_647_merge_other.mat'],'struct_other')
% 
% clear struct_other
% 
% figure('Renderer', 'painters', 'Position', [600 600 275 250])
% hm=heatmap(Median, 'MissingDataColor', [0.8 0.8 0.8]);
% colormap hot
% grid off
% caxis([0 2])
% colorbar off
% %mark the chromosome boundary
% origState = warning('query', 'MATLAB:structOnObject');
% cleanup = onCleanup(@()warning(origState));
% warning('off','MATLAB:structOnObject')
% S = struct(hm); % Undocumented
% ax = S.Axes;    % Undocumented
% clear('cleanup')
% S.XAxis.TickLabelRotation = 45;
% set(gca,'fontsize', 7);
% %set(gca, 'FontName', 'Arial');
% mark = [chrIstop];
% xline(ax, [mark+.5, mark+.5], 'k-');
% yline(ax, [mark+.5, mark+.5], 'k-');
% title(['median spatial distance other cells n=' num2str(n_other)]);
% 
% exportgraphics(gcf, [save_path '/median spatial distance' AgeRange '_othercells.pdf'],'ContentType','vector','BackgroundColor','none');
% 
% figure('Renderer', 'painters', 'Position', [600 600 275 250])
% hm=heatmap(Std, 'MissingDataColor', [0.8 0.8 0.8]);
% colormap parula %abyss
% grid off
% caxis([0 1])
% colorbar off
% %mark the chromosome boundary
% origState = warning('query', 'MATLAB:structOnObject');
% cleanup = onCleanup(@()warning(origState));
% warning('off','MATLAB:structOnObject')
% S = struct(hm); % Undocumented
% ax = S.Axes;    % Undocumented
% clear('cleanup')
% S.XAxis.TickLabelRotation = 45;
% set(gca,'fontsize', 7);
% %set(gca, 'FontName', 'Arial');
% mark = [chrIstop];
% xline(ax, [mark+.5, mark+.5], 'k-');
% yline(ax, [mark+.5, mark+.5], 'k-');
% title(['std other cells'])
% 
% exportgraphics(gcf, [save_path '/std spatial distance' AgeRange '_othercells.pdf'],'ContentType','vector','BackgroundColor','none');
