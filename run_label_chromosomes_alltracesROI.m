clear all
close all

%% load labeled nuclei data
opts = delimitedTextImportOptions("NumVariables", 25);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Volume", "sample_id", "X", "Y", "Z", "ConvexVolume", "EquivDiameter", "Pcell_id", "cell_name", "imageID", "Initial_X", "Initial_Y", "Initial_Z", "Assigned_name", "Aligned_X", "Aligned_Y", "Aligned_Z", "Aligned_X_mod", "Aligned_Y_mod", "Aligned_Z_mod", "min_Procrustes_mod_distance", "min_pattern_Procrustes_mod", "min_Procrustes_distance", "min_pattern_Procrustes", "Ambiguous_state"];
opts.VariableTypes = ["double", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "categorical", "double", "double", "double", "double", "double", "double", "double", "categorical", "double", "categorical", "categorical"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["sample_id", "Assigned_name", "min_pattern_Procrustes_mod", "min_pattern_Procrustes", "Ambiguous_state"], "EmptyFieldRule", "auto");

% Import the data
fulltestdfaligned = readtable("nb04_results/alignment_test_data/df_aligned_XYZ_complete.csv", opts);

clear opts

%% load traces, connect to sample_id, cell_name, and Assigned_name

folder = 'traces_nonp';
Files = dir(['../../traces/' folder '/']); % lists all files

% uncomment 3 lines below if you want to look at only one sample
% Sample = 'Sample7FOV2';
% fulltestdfaligned = fulltestdfaligned(fulltestdfaligned.sample_id == Sample, :);
% Files = dir(['../../traces/' folder '/' Sample '/']); %lists Sample files


mkdir ([folder '_labeled'])
numTraces = length(Files)-2

for p = 3:length(Files)
    File_name = Files(p).name;

    newStr = split(File_name, ["_", "."]);
    TraceSample = [newStr{4:5,1}];
    newStr2=split(File_name, ".");

    load(['../../traces/' folder '/' File_name]) %load territory roi file
    load(['../../Ncoords_1/' TraceSample '.mat']) %load nstats for sample id
    load(['../../raw_images_1/' TraceSample 'ImDAPI.mat']) %load embryo image

    age = height(nstats);

    for j = 1:numel(Trace)
        %for each trace in the territory roi file, can be between 1 and 4
        if numel(Trace) == 1 && iscell(Trace) && length(TAD_id{1})>1 %if there's only one trace in the file the format is different
            roundTrace = round(Trace{1,1});
            for i = 1:height(nstats) %all nuclei in samplefov
                check = ismember(roundTrace,nstats.VoxelList{i,1}, 'rows');%find which nucleus trace belongs to 
                if mean(check) > 0.3
                    TraceNucID = i;
                end
            end
            if exist('TraceNucID','var') == 1
                for q = 1:height(fulltestdfaligned)%find Assigned_name from alignment file using TraceSample and TraceNucID
                    if string(fulltestdfaligned.sample_id(q)) == TraceSample && fulltestdfaligned.cell_name(q) == TraceNucID
                        newStr2{1,1}
                        TraceCellNames(1,1) = string(fulltestdfaligned.Assigned_name(q))%annotate trace with name
                        Procrustes_mod = fulltestdfaligned.min_Procrustes_mod_distance(q);
                        Procrustes = fulltestdfaligned.min_Procrustes_distance(q);
                        TraceCellRecord = q;
                    end
                end
                if p < 50
                    figure
                    imagesc(ImDAPI)
                    colormap gray
                    hold on
                    scatter3(fulltestdfaligned.Initial_X(TraceCellRecord), fulltestdfaligned.Initial_Y(TraceCellRecord), fulltestdfaligned.Initial_Z(TraceCellRecord), 150, 'filled', 'cyan')
                    scatter3(nstats.Centroid(:,1),nstats.Centroid(:,2),nstats.Centroid(:,3), 50, 'filled', 'blue')
                    scatter3(pcellncenter(:,1),pcellncenter(:,2),pcellncenter(:,3), 50, 'filled', 'green')
                    text(pcellncenter(:,1),pcellncenter(:,2)+5,pcellncenter(:,3), 'P-cell')
                    plot3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'magenta', 'LineWidth', 2);
                    text(fulltestdfaligned.Initial_X(TraceCellRecord), fulltestdfaligned.Initial_Y(TraceCellRecord)-5, fulltestdfaligned.Initial_Z(TraceCellRecord),TraceCellNames(1,1))
                    title([newStr{4:7,1} ' label ' TraceCellNames(1,1)])
                    axis equal
                end
            else
                disp([ newStr2{1,1} ' no nucleus found!'])%these are traces outside of nuclei volumes i.e. background or undetectable DAPI
                figure
                imagesc(ImDAPI)
                colormap gray
                hold on
                scatter3(nstats.Centroid(:,1),nstats.Centroid(:,2),nstats.Centroid(:,3), 50, 'filled', 'blue')
                scatter3(pcellncenter(:,1),pcellncenter(:,2),pcellncenter(:,3), 50, 'filled', 'green')
                text(pcellncenter(:,1),pcellncenter(:,2)+5,pcellncenter(:,3), 'P-cell')
                plot3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'magenta', 'LineWidth', 2);
                title([newStr{4:7,1} ' NO LABEL'])
                axis equal
            end

        elseif iscell(TAD_id) && length(TAD_id)>1 %for files with multiple traces in the roi
            %find which nucleus trace belongs to
            roundTrace = round(Trace{1,j});
            if length(roundTrace) >1
                for i = 1:height(nstats) %all nuclei in samplefov
                    check = ismember(roundTrace,nstats.VoxelList{i,1}, 'rows');
                    if mean(check) > 0.3
                        TraceNucID = i;
                    end
                end
                if exist('TraceNucID','var') == 1
                    %find Assigned_name from alignment file using TraceSample and
                    %TraceNucID
                    for q = 1:height(fulltestdfaligned)
                        if string(fulltestdfaligned.sample_id(q)) == TraceSample && fulltestdfaligned.cell_name(q) == TraceNucID
                            newStr2{1,1}
                            TraceCellNames(1,j) = string(fulltestdfaligned.Assigned_name(q))%annotate trace
                            Procrustes_mod = fulltestdfaligned.min_Procrustes_mod_distance(q);
                            Procrustes = fulltestdfaligned.min_Procrustes_distance(q);
                            TraceCellRecord = q;
                        end
                    end
                    if p < 50
                        figure
                        imagesc(ImDAPI)
                        colormap gray
                        hold on
                        scatter3(fulltestdfaligned.Initial_X(TraceCellRecord), fulltestdfaligned.Initial_Y(TraceCellRecord), fulltestdfaligned.Initial_Z(TraceCellRecord), 150, 'filled', 'cyan')
                        scatter3(nstats.Centroid(:,1),nstats.Centroid(:,2),nstats.Centroid(:,3), 50, 'filled', 'blue')
                        scatter3(pcellncenter(:,1),pcellncenter(:,2),pcellncenter(:,3), 50, 'filled', 'green')
                        text(pcellncenter(:,1),pcellncenter(:,2)+5,pcellncenter(:,3), 'P-cell')
                        plot3(Trace{1,j}(:,1), Trace{1,j}(:,2), Trace{1,j}(:,3),'magenta', 'LineWidth', 2);
                        text(fulltestdfaligned.Initial_X(TraceCellRecord), fulltestdfaligned.Initial_Y(TraceCellRecord)-5, fulltestdfaligned.Initial_Z(TraceCellRecord),TraceCellNames(1,j))
                        title([newStr{4:7,1} ' label ' TraceCellNames(1,j)])
                        axis equal
                    end
                else
                    disp([ newStr2{1,1} ' no nucleus found!']) %these are traces outside of nuclei volumes i.e. background or undetectable DAPI
                    figure
                    imagesc(ImDAPI)
                    colormap gray
                    hold on
                    scatter3(nstats.Centroid(:,1),nstats.Centroid(:,2),nstats.Centroid(:,3), 50, 'filled', 'blue')
                    scatter3(pcellncenter(:,1),pcellncenter(:,2),pcellncenter(:,3), 50, 'filled', 'green')
                    text(pcellncenter(:,1),pcellncenter(:,2)+5,pcellncenter(:,3), 'P-cell')
                    plot3(Trace{1,j}(:,1), Trace{1,j}(:,2), Trace{1,j}(:,3),'magenta', 'LineWidth', 2);
                    title([newStr{4:7,1} ' NO LABEL'])
                    axis equal
                end

            else
                disp([ newStr2{1,1} ' less than 1 foci!'])
            end
        end
        %%%%%%%% go to next trace in roi
    end

    if exist('TraceNucID','var') == 1 %save new annotated trace file for all traces in roi
        save([folder '_labeled/' num2str(age) '_' newStr{2,1} '_' newStr{3,1} '_' newStr{4,1} '_' newStr{5,1} '_' newStr{6,1} '_' newStr{7,1} '_labeled.mat'], 'TraceCellNames', 'TraceSample', 'TraceNucID', 'Trace', 'TAD_id','TADsToExclude', 'x1', 'x2', 'y1', 'y3', 'Procrustes_mod', 'Procrustes')
    end
    clearvars -except fulltestdfaligned folder Files
end



