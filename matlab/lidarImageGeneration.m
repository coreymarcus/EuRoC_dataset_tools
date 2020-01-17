% Generates Simulated Lidar measurements from the viewpoint of cam0 in the
% frame of cam0
% Corey Marcus
% UT Austin: Aerospace Engineering

clear
close all
clc

%% Options
LidarFOVHeight = pi/2; %radians
LidarFOVWidth = 3*pi/4; %radians
% LidarFOVHeight = .01; %radians
% LidarFOVWidth = .01; %radians
LidarArrayWidth = 16;
LidarArrayHeight = 16;
pointCloudDownSample = 5; %1 for take every point

%% Main

%load dataset
datasetPath = '~/Documents/EuRoC/V2_01_easy';
addpath('quaternion');
dataset = dataset_load(datasetPath);

%rotations and translations
R_cam2body = dataset.body{1,1}.sensor{1,1}.T_BS(1:3,1:3);
Q_cam2body = dcm2quat(R_cam2body);
P_cam2body = dataset.body{1,1}.sensor{1,1}.T_BS(1:3,4);

%visual-inertial filtered time, position, and quaternions
t_inertial = double(dataset.body{1,1}.sensor{1,5}.data.t)/1E9;
P_inertial2body = dataset.body{1,1}.sensor{1,5}.data.p_RS_R; %position from inertial to body frame in inertial frame
Q_inertial2body = dataset.body{1,1}.sensor{1,5}.data.q_RS; %rotation from inertial to body frame

%point cloud
pointCloud = dataset.body{1,1}.sensor{1,4}.data;
Q = size(pointCloud,2);

%find pose of camera in R frame and camera barrel vector to that frame
N = length(t_inertial);
P_inertial2cam = zeros(3,N);
Q_inertial2cam = zeros(4,N);
Barrel_camInInertial = zeros(3,N);
barrelVect = [0 0 1]';
for ii = 1:N
    P_inertial2cam(:,ii) = quatrotate(quatconj(Q_inertial2body(:,ii)'), P_cam2body')'...
        + P_inertial2body(:,ii);
    Q_inertial2cam(:,ii) = quatmultiply(Q_inertial2body(:,ii)',quatconj(Q_cam2body))';
    Barrel_camInInertial(:,ii) = quatrotate(quatconj(Q_inertial2cam(:,ii)'),barrelVect')';
end

%get camera image times
t_cam = double(dataset.body{1,1}.sensor{1,1}.data.t)/1E9;
L = length(t_cam);

%correlate camera time indicies and inertial time indicies
matchedInertialIdxs = zeros(L,1);
for ii = 1:L
    [~,I] = min(abs(t_inertial - t_cam(ii)));
    matchedInertialIdxs(ii) = I;
end

%build lidar angle array
LidarYawAngles = linspace(-LidarFOVWidth/2,LidarFOVWidth/2,LidarArrayWidth);
LidarPitchAngles = linspace(LidarFOVHeight/2,-LidarFOVHeight/2,LidarArrayHeight); %note intentional sign reversal

%create lidar image arrays
% LidarImages = zeros(LidarArrayHeight,LidarArrayWidth,L);

%build lidar images
tic

LidarImages = lidarImageBuild(LidarArrayHeight,...
    LidarArrayWidth, LidarYawAngles, LidarPitchAngles, 2, Q_inertial2cam,...
    matchedInertialIdxs, P_inertial2cam, pointCloud(1:pointCloudDownSample:end,:));

toc

save('LidarImages','LidarImages')

%% plotting

% figure
% plot3(P_inertial2body(1,:),P_inertial2body(2,:),P_inertial2body(3,:))
% hold on
% plot3(P_inertial2cam(1,:),P_inertial2cam(2,:),P_inertial2cam(3,:))
% scatter3(P_inertial2body(1,1),P_inertial2body(2,1),P_inertial2body(3,1))
% xlabel('x')
% ylabel('y')
% zlabel('z')

figure
plot(P_inertial2body(1,:),P_inertial2body(2,:))
hold on
plot(P_inertial2cam(1,:),P_inertial2cam(2,:))
scatter(P_inertial2body(1,1),P_inertial2body(2,1))
xlabel('x')
ylabel('y')
quiver(P_inertial2cam(1,1:100:end),P_inertial2cam(2,1:100:end),...
    Barrel_camInInertial(1,1:100:end),Barrel_camInInertial(2,1:100:end))
axis equal


figure
plot3(P_inertial2body(1,:),P_inertial2body(2,:),P_inertial2body(3,:))
hold on
plot3(P_inertial2cam(1,:),P_inertial2cam(2,:),P_inertial2cam(3,:))
scatter3(P_inertial2body(1,1),P_inertial2body(2,1),P_inertial2body(3,1))
xlabel('x')
ylabel('y')
zlabel('z')
quiver3(P_inertial2cam(1,1:100:end),P_inertial2cam(2,1:100:end),P_inertial2cam(3,1:100:end),...
    Barrel_camInInertial(1,1:100:end),Barrel_camInInertial(2,1:100:end),Barrel_camInInertial(3,1:100:end))
axis equal

figure
plot3(P_inertial2body(1,:),P_inertial2body(2,:),P_inertial2body(3,:))
hold on
plot3(P_inertial2cam(1,:),P_inertial2cam(2,:),P_inertial2cam(3,:))
scatter3(P_inertial2body(1,1),P_inertial2body(2,1),P_inertial2body(3,1))
scatter3(pointCloud(1:10:end,1),pointCloud(1:10:end,2),pointCloud(1:10:end,3))
xlabel('x')
ylabel('y')
zlabel('z')
quiver3(P_inertial2cam(1,1:100:end),P_inertial2cam(2,1:100:end),P_inertial2cam(3,1:100:end),...
    Barrel_camInInertial(1,1:100:end),Barrel_camInInertial(2,1:100:end),Barrel_camInInertial(3,1:100:end),0)
axis([-5 5 -5 5 0 3])



