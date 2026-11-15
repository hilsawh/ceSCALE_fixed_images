clear all
close all

load('Clustering/AllChromosomes_clusterannotated.mat');

TotalTADNum = 28;
NumClusters = 11;

%% get std values for total and all clusters
AllSTD = [];
STD_total = zeros(TotalTADNum,TotalTADNum);
for i = 1:TotalTADNum
    for j = 1:TotalTADNum
        DisList = [];
        if i ~= j
            for k = 1:length(Chr)
                if Chr(k).r(i) == 1 && Chr(k).r(j) == 1
                    DisList = [DisList ((Chr(k).x(i)-Chr(k).x(j))^2+(Chr(k).y(i)-Chr(k).y(j))^2+(Chr(k).z(i)-Chr(k).z(j))^2)^0.5];
                end
            end
            STD_total(i,j) = std(DisList);
        end
    end
end
STD_total_vector = STD_total(triu(true(size(STD_total)),1));%get the upper triangle of the matrix (no diagonal)
AllSTD = [AllSTD STD_total_vector];

for q = 1:NumClusters
    STD= zeros(TotalTADNum,TotalTADNum);
    for i = 1:TotalTADNum
        for j = 1:TotalTADNum
            DisList = [];
            if i ~= j
                for k = 1:length(Chr)
                    if Chr(k).r(i) == 1 && Chr(k).r(j) == 1
                        if Chr(k).sub0res06sub1res08ID == q
                            DisList = [DisList ((Chr(k).x(i)-Chr(k).x(j))^2+(Chr(k).y(i)-Chr(k).y(j))^2+(Chr(k).z(i)-Chr(k).z(j))^2)^0.5];
                        end
                    end
                end
                STD(i,j) = std(DisList);
            end
        end
    end
    STD_vector = STD(triu(true(size(STD)),1));%get the upper triangle of the matrix (no diagonal)
    AllSTD = [AllSTD STD_vector];
end


%% plot data in boxplot

figure('Renderer', 'painters', 'Position', [50 50 450 250])
boxchart(AllSTD, 'Notch', 'on', 'BoxEdgeColor','k', 'BoxFaceColor', 'w', 'MarkerStyle','.', 'MarkerColor', 'k')
xticklabels([0:11]);
xlabel('cluster')
ylabel('STD of spatial distance')
set(gca,'fontsize', 7);

exportgraphics(gcf, 'STD_total_clusters_box.pdf','ContentType','vector','BackgroundColor','none');
