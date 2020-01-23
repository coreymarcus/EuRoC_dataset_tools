% Convert Generated Lidar Images to a Series of csv files

clear
close all
clc

savedir = 'lidarImages/'; %target directory for lidar images

%load data
load('LidarImagesFinal.mat');

L = size(LidarImages,3);

for ii = 1:L
    filename = strcat(savedir,'lidarImage_',num2str(ii-1,'%05i'),'.csv'); %remember to change padding depending on number of images!
    dlmwrite(filename,LidarImages(:,:,ii),'precision',6);
end