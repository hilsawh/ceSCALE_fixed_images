% This script tests threshold values (LocalMaxThresh) for each probe channel, to be used in
% next program to fit foci in all Hyb rounds
% ------------------------------------------------------------------------------------------------------------------------------------------
% Input: multidimensional matrices (3D image stacks) of foci in a channel
% ------------------------------------------------------------------------------------------------------------------------------------------
% Output: none
% ------------------------------------------------------------------------------------------------------------------------------------------
% Ahilya N. Sawh, PhD
% 06.03.2019
% Version 1.0
% 
% using the Zhuang lab fitMultipleFoci function
%% ------------------------------------------------------------------------------------------------------------------------------------------
clear all
close all

SampleNum = 1;
Hyb = 1;

LocalMaxThresh560 = 70;
LocalMaxThresh647 = 50;
%try otsu or adaptive thresh?

Channel = 560;
FileName = ['sequential/' num2str(Channel) '_' num2str(Hyb) '_' num2str(SampleNum) '.mat'];
load(FileName)

ImageMax = max(ImageStack560,[],3);
figure
imagesc(ImageMax)
colormap gray
axis equal
title('560 max')
figure(100)
imagesc(ImageMax)
colormap gray
axis equal

[Xfit, Yfit, Zfit, Xgof, Ygof, Zgof, Intensity] = fitMultipleFoci(ImageStack560,LocalMaxThresh560);
figure
hold on
plot(Xfit, Yfit, 'ok', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
hold off

figure(100)
hold on
plot(Xfit, Yfit, 'ok', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
for j = 1:length(Xfit)
    text(Xfit(j), Yfit(j), num2str(j));
end
hold off
axis equal
set(gca, 'YDir', 'reverse')
title('560 max foci detected')
saveas(gcf,[num2str(Channel) '_' num2str(Hyb) '_' num2str(SampleNum) '_thresh' num2str(LocalMaxThresh560) '_imagemaxfoci'])

%Pixels Size of Prime95B at 100X (XYZ) (µm):0.11 x 0.11 x 0.20, 1608x1608
%pixels

Zfit = Zfit*0.2; % convert to µm
Xfit = Xfit*110/1000; %µm
Yfit = Yfit*110/1000; %µm
Yfit = 1608*110/1000-Yfit;
    
figure
scatter3(Xfit, Yfit, Zfit, 'ok', 'MarkerFaceColor', 'r');
hold off
axis equal
xlabel('x');
ylabel('y');
zlabel('z');
ylim([0 176.88])
xlim([0 176.88])

%%
Channel = 647;
FileName = ['sequential/' num2str(Channel) '_' num2str(Hyb) '_' num2str(SampleNum) '.mat'];
load(FileName)

ImageMax = max(ImageStack647,[],3);
figure
imagesc(ImageMax)
colormap gray
axis equal
title('647 max')
figure(200)
imagesc(ImageMax)
colormap gray
axis equal


[Xfit, Yfit, Zfit, Xgof, Ygof, Zgof, Intensity] = fitMultipleFoci(ImageStack647,LocalMaxThresh647);
figure
hold on
plot(Xfit, Yfit, 'ok', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
hold off

figure(200)
hold on
plot(Xfit, Yfit, 'ok', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
for j = 1:length(Xfit)
    text(Xfit(j), Yfit(j), num2str(j));
end
hold off
axis equal
set(gca, 'YDir', 'reverse')
title('647 max foci detected')
saveas(gcf,[num2str(Channel) '_' num2str(Hyb) '_' num2str(SampleNum) '_thresh' num2str(LocalMaxThresh647) '_imagemaxfoci'])

%Pixels Size of Prime95B at 100X (XYZ) (µm):0.11 x 0.11 x 0.20, 1608x1608
%pixels

Zfit = Zfit*0.2; % convert to µm
Xfit = Xfit*110/1000; %µm
Yfit = Yfit*110/1000; %µm
Yfit = 1608*110/1000-Yfit;
    
figure
scatter3(Xfit, Yfit, Zfit, 'ok', 'MarkerFaceColor', 'r');
hold off
axis equal
xlabel('x');
ylabel('y');
zlabel('z');
ylim([0 176.88])
xlim([0 176.88])

