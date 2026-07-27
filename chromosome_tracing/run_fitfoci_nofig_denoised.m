function [Xfit, Yfit, Zfit, Xgof, Ygof, Zgof, Intensity] = ANSrun_fitfoci_nofig_denoised(SampleNum)

tic

SampleNum

% assign threshold values for each probe channel
LocalMaxThresh_560 = 10
LocalMaxThresh_647 = 10

TotalNumTADs = 28; %chrI and II together 14 hybs

FileName = ['sequential/560_0_' num2str(SampleNum) '.mat'];
load(['DriftParams' num2str(SampleNum) '.mat']);
load('DeltaZ.mat');
load('tform.mat');
%jetcolors;

load(FileName)



% ImageMax = max(ImageStack,[],3);
% figure(100)
% imagesc(ImageMax)
% colormap gray
% axis equal
% colorbar

% fit  foci

for i = 1:2:TotalNumTADs-1
    j = (i+1)/2;
    
    FileName = ['sequential/647_' num2str(j) '_' num2str(SampleNum) '.mat'];
    load(FileName)
    display('loaded filename')
    
    % warp the 647 image into the 560 channel in xy.
    for m = 1:size(ImageStack647,3)
        ImageStack647(:,:,m) = imtransform(ImageStack647(:,:,m), tform, 'XData', [1 1608], 'Ydata', [1 1608]);
    end
    
    % remove noise
    
    noisyV = ImageStack647;
    
    filteredV = imgaussfilt3(noisyV, 1); %gaussian filter
%     figure
%     subplot(1,2,1)
%     imagesc(max(noisyV,[],3))
%     colormap gray
%     caxis([0 300])
%     subplot(1,2,2)
%     imagesc(max(filteredV,[],3))
%     colormap gray
%     caxis([0 300])
    
    ImageStack647 = filteredV;
    
    display('fitting 647 loci')
    %find and record all foci
    [Xfit{i}, Yfit{i}, Zfit{i}, Xgof{i}, Ygof{i}, Zgof{i}, Intensity{i}] = fitMultipleFoci(ImageStack647,LocalMaxThresh_647);
    Xfit{i} = Xfit{i} - Xdrift(j);
    Yfit{i} = Yfit{i} - Ydrift(j);
    Zfit{i} = Zfit{i} - Zdrift(j);
end

display('fitted 647 loci')

for i = 2:2:TotalNumTADs
    j = i/2;
    
    FileName = ['sequential/560_' num2str(j) '_' num2str(SampleNum) '.mat'];
    load(FileName)
    display('loaded filename')
    
    % remove noise
    
    noisyV = ImageStack560;
    
    filteredV = imgaussfilt3(noisyV, 1); %gaussian filter
%     figure
%     subplot(1,2,1)
%     imagesc(max(noisyV,[],3))
%     colormap gray
%     caxis([0 300])
%     subplot(1,2,2)
%     imagesc(max(filteredV,[],3))
%     colormap gray
%     caxis([0 300])
    
    ImageStack560 = filteredV;
    
    display('fitting 560 loci')
    % find and record all foci
    [Xfit{i}, Yfit{i}, Zfit{i}, Xgof{i}, Ygof{i}, Zgof{i}, Intensity{i}] = fitMultipleFoci(ImageStack560,LocalMaxThresh_560);
    Xfit{i} = Xfit{i} - Xdrift(j);
    Yfit{i} = Yfit{i} - Ydrift(j);
    Zfit{i} = Zfit{i} - Zdrift(j);
end

display('fitted 560 loci')

%% make figures and save data
% figure(200)
% hold on
% for i = 1:TotalNumTADs
%     plot(Xfit{i}, Yfit{i}, 'ok', 'MarkerSize', 5, 'MarkerFaceColor', JetColors{i});
% end
% for i = 1:TotalNumTADs
%     for j = 1:length(Xfit{i})
%         text(Xfit{i}(j), Yfit{i}(j), num2str(i));
%     end
% end
% hold off
% axis equal
% set(gca, 'YDir', 'reverse')

% warp the 647 image into the 560 channel in z.
for i = 1:2:TotalNumTADs-1
    Zfit{i} = Zfit{i}-DeltaZ;
end

% figure(100)
% hold on
% for i = 1:TotalNumTADs
%     scatter3(Xfit{i}, Yfit{i}, Zfit{i}, 'ok');
% end
% hold off
% axis equal
% xlabel('x -px');
% ylabel('y -px');
% zlabel('z -px');
% xlim([0 1608]);
% ylim([0 1608]);
% zlim([0 151]);
% saveas(gcf,['Sample' num2str(SampleNum) '_imagemax_allfoci'])

%Pixels Size of Prime95B at 100X (XYZ) (µm):0.11 x 0.11 x 0.20, 1608x1608
%pixels
for i = 1:TotalNumTADs
    Zfit{i} = Zfit{i}*0.2; %um
    Xfit{i} = Xfit{i}*110/1000; %um
    Yfit{i} = Yfit{i}*110/1000; %um
    Yfit{i} = 1608*110/1000-Yfit{i};
end

% figure(300)
% hold on
% for i = 1:TotalNumTADs
%     scatter3(Xfit{i}, Yfit{i}, Zfit{i}, 'ok', 'MarkerFaceColor', JetColors{i});
% end
% hold off
% axis equal
% xlabel('x -um');
% ylabel('y -um');
% zlabel('z -um');

display('saving result')
save(['result_' num2str(SampleNum) '.mat'], 'Xfit', 'Yfit', 'Zfit','Xgof', 'Ygof', 'Zgof', 'Intensity');

toc
end

