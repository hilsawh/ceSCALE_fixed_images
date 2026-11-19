close all 
clear all
%clearvars -except Chr

load('../AllChromosomes560_647_merge.mat')

res = 'sub0res06_sub1res08';
mkdir (res)
TotalNumClusters = 11; %starting at 1
TotalTADNum = 28;

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 15);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["VarName1", "DNA_trace_snn_res1", "DNA_trace_snn_res11", "DNA_trace_snn_res12", "DNA_trace_snn_res13", "subset_res05_plus_total_identities", "subset_res06_plus_total_identities", "subset_res07_plus_total_identities", "subset_res08_plus_total_identities", "subset_res09_plus_total_identities", "subset_res05_plus_total_identitiesV2", "subset_res06_plus_total_identitiesV2", "subset_res07_plus_total_identitiesV2", "subset_res08_plus_total_identitiesV2", "subset_res09_plus_total_identitiesV2"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
tbl = readtable("/scicore/home/mangos/finatu15/LineageLabeling/Analysis/Traces/chr_matrix_2to150cell/Clustering/Clustering-results_sub0sub1Chr1n2.csv", opts);

% Convert to output type
chromIDs = tbl.VarName1;
clusterIDs_res1 = tbl.DNA_trace_snn_res1;
clusterIDs_res11 = tbl.DNA_trace_snn_res11;
clusterIDs_res12 = tbl.DNA_trace_snn_res12;
clusterIDs_res13 = tbl.DNA_trace_snn_res13;
clusterIDs_sub0res05 = tbl.subset_res05_plus_total_identities;
clusterIDs_sub0res06 = tbl.subset_res06_plus_total_identities;
clusterIDs_sub0res07 = tbl.subset_res07_plus_total_identities;
clusterIDs_sub0res08 = tbl.subset_res08_plus_total_identities;
clusterIDs_sub0res09 = tbl.subset_res09_plus_total_identities;
clusterIDs_sub0res06_sub1res05 = tbl.subset_res05_plus_total_identitiesV2;
clusterIDs_sub0res06_sub1res06 = tbl.subset_res06_plus_total_identitiesV2;
clusterIDs_sub0res06_sub1res07 = tbl.subset_res07_plus_total_identitiesV2;
clusterIDs_sub0res06_sub1res08 = tbl.subset_res08_plus_total_identitiesV2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%this one
clusterIDs_sub0res06_sub1res09 = tbl.subset_res09_plus_total_identitiesV2;

% Clear temporary variables
clear opts tbl

%% annotate total Chr structure with cluster id 
names =  [0, 0.1, 0.2, 1, 1.1, 1.2, 2, 3, 4, 5, 6]; %will be converted to whole integers for simplicity

for q = 1:TotalNumClusters    
    for i = 1:length(clusterIDs_sub0res06_sub1res08) %choose resolution here and next line
        if clusterIDs_sub0res06_sub1res08(i,1) == names(1,q) && iscell(clusterIDs_sub0res06_sub1res08(i,1)) == 0 
            Chr(chromIDs(i,1)).sub0res06sub1res08ID = q; %clusterIDs_sub0res06_sub4res08(i,1);
        end
    end
end
save('AllChromosomes_clusterannotated.mat', 'Chr')


%% create new Chr structure for each cluster and plot

for q = 1:TotalNumClusters
    n = 0;
    for i = 1:length(clusterIDs_sub0res06_sub1res08) %choose resolution here and next line
        if clusterIDs_sub0res06_sub1res08(i,1) == names(1,q) && iscell(clusterIDs_sub0res06_sub1res08(i,1)) == 0 
            n = n+1;
            ChrClus(n).x = Chr(chromIDs(i,1)).x;
            ChrClus(n).y = Chr(chromIDs(i,1)).y;
            ChrClus(n).z = Chr(chromIDs(i,1)).z;
            ChrClus(n).r = Chr(chromIDs(i,1)).r;
            ChrClus(n).Age = Chr(chromIDs(i,1)).Age;
            ChrClus(n).TraceNo = Chr(chromIDs(i,1)).TraceNo;
            ChrClus(n).SampleID = Chr(chromIDs(i,1)).SampleID;
            ChrClus(n).Dataset = Chr(chromIDs(i,1)).Dataset;
            ChrClus(n).CellID = Chr(chromIDs(i,1)).CellID;
            ChrClus(n).NucID = Chr(chromIDs(i,1)).NucID;
            ChrClus(n).Procrustes = Chr(chromIDs(i,1)).Procrustes;
            ChrClus(n).Procrustes_mod = Chr(chromIDs(i,1)).Procrustes_mod;
        end
    end
    
    save ([res '/ChrClus' num2str(q) '.mat'],'ChrClus');
    
    Mean = zeros(TotalTADNum,TotalTADNum);
    Std = zeros(TotalTADNum,TotalTADNum);
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
            Std(i,j) = std(DisList);
            SEM(i,j) = std(DisList)/(length(DisList))^0.5;
            NofData(i,j) = length(DisList);
            DisListAll{i,j} = DisList;
            Median(i,j) = median(DisList);
        end
    end

    save ([res '/DisListClus' num2str(q) '.mat'],'DisListAll');


    % make figures
    
    figure('Renderer', 'painters', 'Position', [550 550 250 200])
    h = heatmap(Median, 'MissingDataColor',[0.8 0.8 0.8]);
    grid off
    s = struct(h);
    s.XAxis.TickLabelRotation = 45;
    set(gca,'fontsize', 7);
    set(gca, 'FontName', 'Arial');
    colorbar
    title(['Median Cluster ' num2str(q) ' res' num2str(res) ' n = ' num2str(n)]);
    colormap hot
    caxis([0 2])
    print([res '/Median_cluster' num2str(q)  ' n = ' num2str(n)],'-dsvg');
    saveas(gcf,[res '/Median_cluster' num2str(q) ' n = ' num2str(n)]);
    exportgraphics(gcf, [res '/Median_cluster' num2str(q)  ' n = ' num2str(n) '.pdf'],'ContentType','vector','BackgroundColor','none');
    
    figure('Renderer', 'painters', 'Position', [550 550 250 200])
    imagesc(Std)
    title(['Std of pairwise distance (um) cluster' num2str(q) res ' n = ' num2str(n)]);
    colorbar
    caxis([0 1.2])
    axis square
    print([res '/Std_cluster' num2str(q)  ' n = ' num2str(n)],'-dsvg');
    saveas(gcf,[res '/Std_cluster' num2str(q) ' n = ' num2str(n)]);
    
    lengths = [];
    for i = 1:length(ChrClus)
        newlength = sum(sum(ChrClus(i).r));
        lengths = [lengths newlength];
    end
    
    figure('Renderer', 'painters', 'Position', [550 550 250 200])
    histogram(lengths)
    xlabel('length of trace')
    ylabel('number of traces in cluster')
    axis square
    xlim([0 28.5])
    title(['Lengths cluster' num2str(q) res ' n = ' num2str(n)]);
    print([res '/Lengths_cluster' num2str(q)  ' n = ' num2str(n)], '-dsvg');
    saveas(gcf,[res '/Lengths_cluster' num2str(q) ' n = ' num2str(n)]);
    
    figure('Renderer', 'painters', 'Position', [550 550 250 200])
    histogram([ChrClus.Age])
    xlabel('age of trace')
    ylabel('number of traces in cluster')
    axis square
    %xlim([0 40.5])
    title(['Ages cluster' num2str(q) res ' n = ' num2str(n)]);
    print([res '/Ages_cluster' num2str(q) ' n = ' num2str(n)], '-dsvg');
    saveas(gcf,[res '/Ages_cluster' num2str(q) ' n = ' num2str(n)]);
    
    clear ChrClus
end
