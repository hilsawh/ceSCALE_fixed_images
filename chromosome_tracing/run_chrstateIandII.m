clear all
close all

save_dir = 'figs';
%% Import data from text file

opts = delimitedTextImportOptions("NumVariables", 4);

opts.DataLines = [1, Inf];
opts.Delimiter = ";";

opts.VariableNames = ["chromosome", "start", "stop", "stateID"];
opts.VariableTypes = ["categorical", "double", "double", "double"];

opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

opts = setvaropts(opts, "chromosome", "EmptyFieldRule", "auto");

tbl = readtable("Evans_EEautosome_states.csv", opts);

chromosome = tbl.chromosome;
start = tbl.start;
stop = tbl.stop;
stateID = tbl.stateID;

clear opts

%%
numstates = 20;
n = 28; %number regions

allstates = zeros(numstates, n);

for state = 1:numstates
%% calculate enrichment
    chrom_both = 'IandII';

    fid = fopen(['chr' chrom_both 'TadDomainStartsAndEnds.txt'], 'r');
    tline = fgetl(fid);
    i = 0;
    while ischar(tline)
        i = i+1;
        C = strsplit(tline);
        DomainStart(i) = str2num(C{2});
        DomainEnd(i) = str2num(C{3});
        DomainCenter(i) = (DomainStart(i)+DomainEnd(i))/2;
        DomainSize(i) = DomainEnd(i)-DomainStart(i);
        tline = fgetl(fid);
    end
    fclose(fid);

    %count number of hits per region chrI
    HperRegion = zeros(1,n);
    chrom = 'I';

    start_filt = start(chromosome==chrom,:);
    start_filt_names = stateID(chromosome==chrom,:);

    for i = 1:14

        for j = 1:length(start_filt)
            if start_filt(j,1) >= DomainStart(1,i) && start_filt(j,1) <= DomainEnd(1,i) && start_filt_names(j,1) == state
                HperRegion(1,i) = (HperRegion(1,i) + 1);
            end
        end

    end

    %count number of hits per region chrII
    chrom = 'II';

    start_filt = start(chromosome==chrom,:);
    start_filt_names = stateID(chromosome==chrom,:);

    for i = 15:28

        for j = 1:length(start_filt)
            if start_filt(j,1) >= DomainStart(1,i) && start_filt(j,1) <= DomainEnd(1,i) && start_filt_names(j,1) == state
                HperRegion(1,i) = (HperRegion(1,i) + 1);
            end
        end

    end

    
    %NormHperRegion=(HperRegion/sum(HperRegion))*100; %normalize only to regions of interest
    %HperRegion=(HperRegion/length(start_filt))*100; %normalize to whole chr data

    % figure
    % bar(HperRegion, 'LineWidth',1);
    % xlabel('region #')
    % ylabel(['# state ' num2str(state)])
    % title(['# state ' num2str(state) ' chr' chrom_both])
    % xlim([0.5 n+0.5])

   
    allstates(state,:) = HperRegion;


    % figure('Renderer', 'painters', 'Position', [250 250 600 100])
    % plot(NormHperRegion, 'LineWidth',2, 'Color', 'k');
    % %xlabel('region')
    % ylabel(['state ' num2str(state)])
    % %title(['normalized state ' num2str(state) ' chr' chrom_both])
    % xlim([0.5 n+0.5])
    % set(gca,'fontsize', 7);
    % set(gca, 'XTickLabel', [])
    % box off
    % 
    % exportgraphics(gcf, [save_dir '/norm state ' num2str(state) ' chr' chrom_both '.pdf'],'ContentType','vector','BackgroundColor','none');

end
writematrix(allstates, [save_dir '/allstates.csv'])
Norm_allstates = normr(allstates);%row normalize

%% plot heatmap of all states

figure('Renderer', 'painters', 'Position', [250 250 600 300])
h = heatmap(Norm_allstates);
colormap(flipud(pink))
%or load custom colormap
load([save_dir '/grayredcmap.mat'])
colormap(CustomColormap)
grid off
%s = struct(h);
%s.XAxis.TickLabelRotation = 45;
set(gca,'fontsize', 7);
xlabel('region')
ylabel('state')
%colorbar off
exportgraphics(gcf, [save_dir '/norm allstates chr' chrom_both '.pdf'],'ContentType','vector','BackgroundColor','none');

%% plot euclidean clustergram https://www.mathworks.com/help/bioinfo/ug/working-with-the-clustergram-function.html
pats = {'1- promoter' '2- 5proximal and gene body' '3- txn elongation I exon' '4- txn elongation II exon and intron' ...
    '5- txn elongation III exon and gene end' '6- txn elongation IV low exp and repeats' '7- txn elongation V introns and repeats'...
    '8- enhancer I intronic' '9- enhancer II intergenic' '10- enhancer III weak' '11- border' '12- repeats, low exp introns, ATrich'...
    '13- retrotransposons, pseudogenes, H3K9me3, H3K27me3' '14- mixed, tissue specific' '15- repeats, RNA pseudogenes, H3K9me2' ...
    '16- intergenic, silent genes, piRNAs and repeats' '17- Pc/H3K27me3 I:low expr/silent and pseudogenes' '18- Pc/H3K27me3 II: low expr/silent' ...
    '19- Pc/H3K27me3 II: low expr/silent, gene body' '20- H3K9me3 and H3K27me3: silent genes and pseudogenes' };
cg_s = clustergram(Norm_allstates, 'RowLabels', pats,...
                               'Cluster', 'Column',...
                               'ImputeFun', @knnimpute)

cg_s.ColumnPDist = 'euclidean';
cg_s.Colormap = redbluecmap;
cg_s.DisplayRange = 200; % max(max(allstates));
set(0,'ShowHiddenHandles','on')
colorbar west
% Get all handles from root
allhnds = get(0,'Children');
% Find heatmap axis and change the font size
h = findall(allhnds, 'Tag', 'HeatMapAxes');
set(h, 'FontSize', 7)

set(gcf, 'Position', [250 250 600 300])

exportgraphics(gcf, [save_dir '/norm allstates chr' chrom_both '_clustergram.pdf'],'ContentType','vector','BackgroundColor','none');

%% plot cosine clustergram https://www.mathworks.com/help/bioinfo/ug/working-with-the-clustergram-function.html

cg_s = clustergram(Norm_allstates, 'RowLabels', pats,...
                               'Cluster', 'Column',...
                               'ImputeFun', @knnimpute, 'ColumnPDist', 'cosine')

cg_s.ColumnPDist = 'cosine';
cg_s.Colormap = redbluecmap;
cg_s.DisplayRange = 200; % max(max(allstates));
set(0,'ShowHiddenHandles','on')
colorbar west
% Get all handles from root
allhnds = get(0,'Children');
% Find heatmap axis and change the font size
h = findall(allhnds, 'Tag', 'HeatMapAxes');
set(h, 'FontSize', 7)

set(gcf, 'Position', [250 250 600 300])

exportgraphics(gcf, [save_dir '/norm allstates cosine chr' chrom_both '.pdf'],'ContentType','vector','BackgroundColor','none');

%% correlate states with mean metaloci signal for AB lineage

clearvars -except allstates Norm_allstates save_dir numstates n
pat = 'AB';
type = 'hl';

mean_percent = load(['/scicore/home/mangos/GROUP/Kostas_Ahilya_shared_folder/figures_metaloci_percentages/figures_zygotic/mean_percents_' pat '_' type '.csv' ]);

Pearson = zeros(numstates, 1);
for i = 1:numstates
    R = corr(mean_percent', Norm_allstates(i,:)');
    Pearson(i,1) = R;
end

figure('Renderer', 'painters', 'Position', [250 250 400 100])
bar(Pearson, 'k')
box off
ylim([-1 1])
xlim([0.5 20.5])
set(gca, 'FontSize', 7)
ylabel(['R ' type])
xlabel('state')

exportgraphics(gcf, [save_dir '/pearson_state_' pat '_' type '.pdf'],'ContentType','vector','BackgroundColor','none');
