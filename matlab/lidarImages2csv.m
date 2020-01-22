% Convert Generated Lidar Images to a Series of csv files

clear
close all
clc

savedir = 'lidarImages/'; %target directory for lidar images

%load data
load('LidarImages.mat');

L = size(LidarImages,3);

for ii = 1:L
    filename = strcat(savedir,'lidarImage_',string(ii-1),'.csv');
    dlmwrite(filename,LidarImages(:,:,ii),'precision',6);
end