tic;
clear all
close all

Folder{1} = '200626_traces';
Folder{2} = '200702_traces';
Folder{3} = '201123_traces';
Folder{4} = '201124_traces';
Traces_folder = {'traces_nonp_labeled', 'traces_p_labeled'};


AgeStart = 2; %change these to specify embryo age
AgeStop = 150;
AgeRange = ([ num2str(AgeStart) 'to' num2str(AgeStop) 'cell']);

SegChannel = 647;
save_path = ['chr_matrix_' AgeRange];
mkdir(save_path)
RegionStart = 1;
TotalTADNum = 28;
n = 0;
flag = 0;

for ii = 1:length(Folder)
    for i = 1:length(Traces_folder)
        [Folder{ii} '/' Traces_folder{i}]

        Files = dir([Folder{ii} '/' Traces_folder{i}]); % lists all files
        for p = 1:length(Files)
            File_name = Files(p).name;
            newStr = split(File_name, ["roi", "_"]);


            q = str2num(strtok(File_name,'_'));
            if isempty(q)==0 && AgeStart <= q && q <= AgeStop && sum(newStr{7,1} == num2str(SegChannel))>0 % select ages of interest
                load([Folder{ii} '/' Traces_folder{i} '/' File_name]);
                TraceName = [Folder{ii} '/' Traces_folder{i} '/' File_name];
                flag = 1;
            end

            if flag == 1
                if length(Trace) ==4
                    for g = 1:4
                        if length(TAD_id{g})>1 %&& length(TAD_id{1,g}) == length(unique(TAD_id{1,g})) %skip trace if it has multiple foci per TAD/region (remove ambiguous traces)
                            n = n+1;
                            Chr(n).x = zeros(TotalTADNum,1);
                            Chr(n).y = zeros(TotalTADNum,1);
                            Chr(n).z = zeros(TotalTADNum,1);
                            Chr(n).r = zeros(TotalTADNum,1);

                            for k = 1: length(TAD_id{1,g})
                                Chr(n).x(TAD_id{1,g}(k,1)) = Trace{1,g}(k,1)*110/1000; % convert to um from pixel units
                                Chr(n).y(TAD_id{1,g}(k,1)) = 1608*110/1000-Trace{1,g}(k,2)*110/1000;
                                Chr(n).z(TAD_id{1,g}(k,1)) = Trace{1,g}(k,3)*0.2;
                                Chr(n).r(TAD_id{1,g}(k,1)) = 1;
                                Trace{1,g}(k,1) = Trace{1,g}(k,1)*110/1000;
                                Trace{1,g}(k,2) = Trace{1,g}(k,2)*110/1000;
                                Trace{1,g}(k,3) = Trace{1,g}(k,3)*0.2;
                            end
                            
                            Chr(n).Age = q;
                            Chr(n).TraceNo = g;
                            Chr(n).TracesInTerritory = 4;
                            Chr(n).Procrustes = Procrustes;
                            Chr(n).Procrustes_mod = Procrustes_mod;
                            Chr(n).TraceName = TraceName;
                            Chr(n).SampleID  = TraceSample;
                            Chr(n).Dataset = Folder{ii};
                            Chr(n).CellID  = TraceCellNames{1,g};
                            Chr(n).NucID  = TraceNucID;
                        end
                    end
                end
                if length(Trace) ==3
                    for g = 1:3
                        if length(TAD_id{g})>1 %&& length(TAD_id{1,g}) == length(unique(TAD_id{1,g}))
                            n = n+1;
                            Chr(n).x = zeros(TotalTADNum,1);
                            Chr(n).y = zeros(TotalTADNum,1);
                            Chr(n).z = zeros(TotalTADNum,1);
                            Chr(n).r = zeros(TotalTADNum,1);

                            for k = 1: length(TAD_id{1,g})
                                Chr(n).x(TAD_id{1,g}(k,1)) = Trace{1,g}(k,1)*110/1000; % convert to um from pixel units
                                Chr(n).y(TAD_id{1,g}(k,1)) = 1608*110/1000-Trace{1,g}(k,2)*110/1000;
                                Chr(n).z(TAD_id{1,g}(k,1)) = Trace{1,g}(k,3)*0.2;
                                Chr(n).r(TAD_id{1,g}(k,1)) = 1;
                                Trace{1,g}(k,1) = Trace{1,g}(k,1)*110/1000;
                                Trace{1,g}(k,2) = Trace{1,g}(k,2)*110/1000;
                                Trace{1,g}(k,3) = Trace{1,g}(k,3)*0.2;
                            end
                            Chr(n).Age = q;
                            Chr(n).TraceNo = g;
                            Chr(n).TracesInTerritory = 3;
                            Chr(n).Procrustes = Procrustes;
                            Chr(n).Procrustes_mod = Procrustes_mod;
                            Chr(n).TraceName = TraceName;
                            Chr(n).SampleID  = TraceSample;
                            Chr(n).Dataset = Folder{ii};
                            Chr(n).CellID  = TraceCellNames{1,g};
                            Chr(n).NucID  = TraceNucID;
                        end
                    end
                end
                if length(Trace) ==2
                    for g = 1:2
                        if length(TAD_id{g})>1 %&& length(TAD_id{1,g}) == length(unique(TAD_id{1,g}))
                            n = n+1;
                            Chr(n).x = zeros(TotalTADNum,1);
                            Chr(n).y = zeros(TotalTADNum,1);
                            Chr(n).z = zeros(TotalTADNum,1);
                            Chr(n).r = zeros(TotalTADNum,1);

                            for k = 1: length(TAD_id{1,g})
                                Chr(n).x(TAD_id{1,g}(k,1)) = Trace{1,g}(k,1)*110/1000; % convert to um from pixel units
                                Chr(n).y(TAD_id{1,g}(k,1)) = 1608*110/1000-Trace{1,g}(k,2)*110/1000;
                                Chr(n).z(TAD_id{1,g}(k,1)) = Trace{1,g}(k,3)*0.2;
                                Chr(n).r(TAD_id{1,g}(k,1)) = 1;
                                Trace{1,g}(k,1) = Trace{1,g}(k,1)*110/1000;
                                Trace{1,g}(k,2) = Trace{1,g}(k,2)*110/1000;
                                Trace{1,g}(k,3) = Trace{1,g}(k,3)*0.2;
                            end
                            
                            Chr(n).Age = q;
                            Chr(n).TraceNo = g;
                            Chr(n).TracesInTerritory = 2;
                            Chr(n).Procrustes = Procrustes;
                            Chr(n).Procrustes_mod = Procrustes_mod;
                            Chr(n).TraceName = TraceName;
                            Chr(n).SampleID  = TraceSample;
                            Chr(n).Dataset = Folder{ii};
                            Chr(n).CellID  = TraceCellNames{1,g};
                            Chr(n).NucID  = TraceNucID;
                        end
                    end
                end
                if length(Trace) ==1 && iscell(Trace) && length(TAD_id{1})>1 %&& length(TAD_id{1,1}) == length(unique(TAD_id{1,1}))
                    n = n+1;
                    Chr(n).x = zeros(TotalTADNum,1);
                    Chr(n).y = zeros(TotalTADNum,1);
                    Chr(n).z = zeros(TotalTADNum,1);
                    Chr(n).r = zeros(TotalTADNum,1);

                    for k = 1: length(TAD_id{1,1})
                        Chr(n).x(TAD_id{1,1}(k,1)) = Trace{1,1}(k,1)*110/1000; % convert to um from pixel units
                        Chr(n).y(TAD_id{1,1}(k,1)) = 1608*110/1000-Trace{1,1}(k,2)*110/1000;
                        Chr(n).z(TAD_id{1,1}(k,1)) = Trace{1,1}(k,3)*0.2;
                        Chr(n).r(TAD_id{1,1}(k,1)) = 1;
                        Trace{1,1}(k,1) = Trace{1,1}(k,1)*110/1000;
                        Trace{1,1}(k,2) = Trace{1,1}(k,2)*110/1000;
                        Trace{1,1}(k,3) = Trace{1,1}(k,3)*0.2;
                    end
                    
                    Chr(n).Age = q;
                    Chr(n).TraceNo = 1;
                    Chr(n).TracesInTerritory = 1;
                    Chr(n).Procrustes = Procrustes;
                    Chr(n).Procrustes_mod = Procrustes_mod;
                    Chr(n).TraceName = TraceName;
                    Chr(n).SampleID  = TraceSample;
                    Chr(n).Dataset = Folder{ii};
                    Chr(n).CellID  = TraceCellNames;
                    Chr(n).NucID  = TraceNucID;
                end

                flag = 0;
            end
        end
    end
end

%%
Mean = zeros(TotalTADNum,TotalTADNum);
Std = zeros(TotalTADNum,TotalTADNum);
Median = zeros(TotalTADNum,TotalTADNum);


for i = 1:TotalTADNum
    for j = 1:TotalTADNum
        if i ~=j
            DisList = [];
            for k = 1:length(Chr)
                if Chr(k).r(i) == 1 && Chr(k).r(j) == 1
                    if ((Chr(k).x(i)-Chr(k).x(j))^2+(Chr(k).y(i)-Chr(k).y(j))^2+(Chr(k).z(i)-Chr(k).z(j))^2)^0.5 <=6
                        DisList = [DisList ((Chr(k).x(i)-Chr(k).x(j))^2+(Chr(k).y(i)-Chr(k).y(j))^2+(Chr(k).z(i)-Chr(k).z(j))^2)^0.5];
                    end
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


%% make figures

figure('Renderer', 'painters', 'Position', [150 150 210 200])
h = heatmap(Mean(RegionStart:TotalTADNum,RegionStart:TotalTADNum));
originalUnits = h.Units;
h.Units = 'centimeters';
sz = size(h.ColorData);
originalPos = h.Position;
h.Position(3:4) = min(h.Position(3:4))*[1,1];
grid off
colormap hot
caxis([0 2])
% h.XData = [round(DomainCenter(1,RegionStart:end),1)];
% h.YData = [round(DomainCenter(1,RegionStart:end),1)];
title('mean spatial distance')
set(gca, 'FontName', 'Arial', 'FontSize', 7)
%saveas(gcf, 'mean');
print([save_path '/' 'mean_' num2str(SegChannel) ], '-dsvg');
exportgraphics(gcf, [save_path '/mean spatial distance' AgeRange '_' SegChannel '.pdf'],'ContentType','vector','BackgroundColor','none');

figure('Renderer', 'painters', 'Position', [150 150 210 200])
h = heatmap(Median(RegionStart:TotalTADNum,RegionStart:TotalTADNum));
originalUnits = h.Units;
h.Units = 'centimeters';
sz = size(h.ColorData);
originalPos = h.Position;
h.Position(3:4) = min(h.Position(3:4))*[1,1];
grid off
colormap hot
caxis([0 2])
% h.XData = [round(DomainCenter(1,RegionStart:end),1)];
% h.YData = [round(DomainCenter(1,RegionStart:end),1)];
title('median spatial distance')
set(gca, 'FontName', 'Arial', 'FontSize', 7)
%saveas(gcf, 'median');
print([save_path '/''median' num2str(SegChannel) ], '-dsvg');
exportgraphics(gcf, [save_path '/median spatial distance' AgeRange '_' SegChannel '.pdf'],'ContentType','vector','BackgroundColor','none');

figure
imagesc(SEM)
colorbar
title(['SEM of pairwise distance - ' AgeRange]);
xticks([RegionStart:TotalTADNum])
% xticklabels({round(DomainCenter(1,RegionStart:end),1)})
yticks([RegionStart:TotalTADNum])
% yticklabels({round(DomainCenter(1,RegionStart:end),1)})
xtickangle(45)
%PlotProp
axis square
%saveas(gcf, 'sem');
print([save_path '/''SEM' num2str(SegChannel) ],'-dsvg')
exportgraphics(gcf, [save_path '/SEM' AgeRange '_' SegChannel '.pdf'],'ContentType','vector','BackgroundColor','none');

figure
imagesc(Std)
colorbar
title(['Std of pairwise distance - ' AgeRange]);
xticks([RegionStart:TotalTADNum])
% xticklabels({round(DomainCenter(1,RegionStart:end),1)})
yticks([RegionStart:TotalTADNum])
% yticklabels({round(DomainCenter(1,RegionStart:end),1)})
xtickangle(45)
%PlotProp
axis square
%saveas(gcf, 'std');
print([save_path '/''STD' num2str(SegChannel) ],'-dsvg')
exportgraphics(gcf, [save_path '/STD' AgeRange  '_' SegChannel '.pdf'],'ContentType','vector','BackgroundColor','none');


figure
imagesc(NofData)
colorbar
title(['Number of measurements - ' AgeRange]);
xticks([RegionStart:TotalTADNum])
% xticklabels({round(DomainCenter(1,RegionStart:end),1)})
yticks([RegionStart:TotalTADNum])
% yticklabels({round(DomainCenter(1,RegionStart:end),1)})
xtickangle(45)
%PlotProp
axis square
%saveas(gcf, 'nofdata');
print([save_path '/''NofDATA' num2str(SegChannel) ],'-dsvg')
exportgraphics(gcf, [save_path '/NofData' AgeRange '_' SegChannel '.pdf'],'ContentType','vector','BackgroundColor','none');


% figure
% histogram([Chr.TracesInTerritory])
% title('wild-type segmentation results')
% ylabel('number of traces')
% xlabel('total traces in territory volume')
%PlotProp

TracesInTerritory = [Chr.TracesInTerritory].';

TracesPerSegment = [];
TracesPerSegment(1,1) = sum(TracesInTerritory(:) == 1)/length(Chr)*100;
TracesPerSegment(1,2) = sum(TracesInTerritory(:) == 2)/length(Chr)*100;
TracesPerSegment(1,3) = sum(TracesInTerritory(:) == 3)/length(Chr)*100;
TracesPerSegment(1,4) = sum(TracesInTerritory(:) == 4)/length(Chr)*100;

figure
bar(TracesPerSegment)
title('wild-type segmentation results')
ylabel('% of total traces')
xlabel('traces within territory volume')
ylim([0 100])
%PlotProp
print([save_path '/segmentation' num2str(SegChannel) ],'-dsvg')
exportgraphics(gcf, [save_path '/segmentation' AgeRange  '_' SegChannel '.pdf'],'ContentType','vector','BackgroundColor','none');


save([save_path '/MeanPairDistanceMatrix' num2str(SegChannel) '.mat'], 'Mean');
save([save_path '/MedianPairDistanceMatrix' num2str(SegChannel) '.mat'], 'Median');
save([save_path '/StdPairDistanceMatrix' num2str(SegChannel) '.mat'], 'Std');
save([save_path '/SEMPairDistanceMatrix' num2str(SegChannel) '.mat'], 'SEM');
save([save_path '/NofDataPairDistanceMatrix' num2str(SegChannel) '.mat'], 'NofData');
save([save_path '/AllChromosomes' num2str(SegChannel) '.mat'], 'Chr');
save([save_path '/DisListAll' num2str(SegChannel) '.mat'], 'DisListAll');

toc