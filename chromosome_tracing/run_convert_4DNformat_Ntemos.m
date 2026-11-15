clear all
close all

save_path = 'Ntemos2025/';

%% import traces in matlab struct format and coordinates

load([save_path 'AllChromosomes_clusterannotated.mat']);

opts = spreadsheetImportOptions("NumVariables", 3);
opts.Sheet = "probe locations";
opts.DataRange = "A2:C29";
opts.VariableNames = ["chromosome", "start", "stop"];
opts.VariableTypes = ["categorical", "double", "double"];
opts = setvaropts(opts, "chromosome", "EmptyFieldRule", "auto");
coords = readtable([save_path '/Supplementary Data 9.xlsx'], opts, "UseExcel", false);
clear opts
%% run conversion for each dataset at a time
dataset = Chr;
name = 'chrIchrII';


T = table;
n = 0;
for i = 1:length(dataset) %total traces
    for j = 1:length(dataset(i).r) %number of regions probed in experiment
        if dataset(i).r(j,1) == 1 %if region was detected
            
            n = n+1;
            T.Spot_ID(n) = n;
            T.Trace_ID(n) = i;
            T.X(n) = dataset(i).x(j,1);
            T.Y(n) = dataset(i).y(j,1);
            T.Z(n) = dataset(i).z(j,1);         
            T.Chrom(n) = coords.chromosome(j);
            T.Chrom_Start(n) = coords.start(j);
            T.Chrom_End(n) = coords.stop(j);

        end 
    end
end

writetable(T, [save_path name '4DN.csv'],'Delimiter',',','QuoteStrings','none')
