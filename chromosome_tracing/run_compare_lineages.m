%% compare average matrices of chromosome conformation in different lineages
clear all 
close all

pats = {'ABa' 'ABp' 'E' 'MS' 'C' 'D' 'P'};
numpats = size(pats, 2);

TotalTADNumBoth = 28;
chrIstart = 1;
chrIstop = 14;
chrIIstart = 15;
chrIIstop = 28;

count = 0;
figure ('Renderer', 'painters', 'Position', [10 10 1200 1200])

for a = 1:numpats
    for b = 1:numpats
        %load dataset1, get median
        load(['AllChromosomes560_647_merge_'  pats{1,a} '.mat'])
        mx1 = zeros(TotalTADNumBoth,TotalTADNumBoth);
        for i = chrIstart:chrIIstop
            for j = chrIstart:chrIIstop
                if i ~= j
                    DisList = [];
                    for k = 1:length(struct)
                        if struct(k).r(i) == 1 && struct(k).r(j) == 1 %&& struct(k).TracesInTerritory == 1    %pick only well-separated chromosomes
                            DisList = [DisList ((struct(k).x(i)-struct(k).x(j))^2+(struct(k).y(i)-struct(k).y(j))^2+(struct(k).z(i)-struct(k).z(j))^2)^0.5];
                        end
                    end
                    mx1(i,j) = median(DisList);
                end
            end
        end
        %load dataset2, get median
        load(['AllChromosomes560_647_merge_'  pats{1,b} '.mat'])
        mx2 = zeros(TotalTADNumBoth,TotalTADNumBoth);
        for i = chrIstart:chrIIstop
            for j = chrIstart:chrIIstop
                if i ~= j
                    DisList = [];
                    for k = 1:length(struct)
                        if struct(k).r(i) == 1 && struct(k).r(j) == 1 %&& struct(k).TracesInTerritory == 1    %pick only well-separated chromosomes
                            DisList = [DisList ((struct(k).x(i)-struct(k).x(j))^2+(struct(k).y(i)-struct(k).y(j))^2+(struct(k).z(i)-struct(k).z(j))^2)^0.5];
                        end
                    end
                    mx2(i,j) = median(DisList);
                end
            end
        end
        %create difference matrix
        diffmx = mx1 -mx2;

        %plot heatmap
        count = count + 1;
        subplot(numpats, numpats, count)
        hm = heatmap(diffmx);
        grid off
        caxis([-1 1])
        colormap redbluecmap
        colorbar off
        % %mark the chromosome boundary
        % origState = warning('query', 'MATLAB:structOnObject');
        % cleanup = onCleanup(@()warning(origState));
        % warning('off','MATLAB:structOnObject')
        % S = struct(hm); % Undocumented
        % ax = S.Axes;    % Undocumented
        % clear('cleanup')
        % %S.XAxis.TickLabelRotation = 45;
        % set(gca,'fontsize', 7);
        % mark = 14;
        % xline(ax, [mark+.5, mark+.5], 'k-');
        % yline(ax, [mark+.5, mark+.5], 'k-');
        title([pats{1,a} '-' pats{1,b}]);

        clear mx1 mx2 diffmx


    end
end

exportgraphics(gcf, ['Lineage2LineageDiffs.pdf'],'ContentType','vector','BackgroundColor','none');
saveas(figure(1),['Lineage2LineageDiffs'],'fig');
