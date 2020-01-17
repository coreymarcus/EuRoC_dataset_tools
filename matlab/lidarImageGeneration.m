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
LidarArrayWidth = 64;
LidarArrayHeight = 64;

%% Main

%load dataset
datasetPath = '~/Documents/EuRoC/V2_01_easy';
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
LidarImages = zeros(LidarArrayHeight,LidarArrayWidth,L);

%build lidar images
tic
for ii = 1:1 %iteration on images
    parfor jj = 1:LidarArrayHeight
        for kk = 1:LidarArrayWidth
%             tic
            %create quaternion corresponding to this angle
            Q_cam2lidarray = angle2quat(LidarYawAngles(kk),LidarPitchAngles(jj),0,'ZXY');
            
            %generate a lidar ray in the inertial frame
            Q_inertial2lidarray = quatmultiply(Q_inertial2cam(:,matchedInertialIdxs(ii))',...
                Q_cam2lidarray);
            
            %transform point cloud to lidar frame
            pointCloudTrans = quatrotate(Q_inertial2lidarray, ...
                pointCloud - P_inertial2cam(:,matchedInertialIdxs(ii))');
%             toc
            
            %find valid indicies in the point cloud trans
%             tic
            validIdxs = (pointCloudTrans(:,3) > .1) & (abs(pointCloudTrans(:,1)) < 0.01) ...
                & (abs(pointCloudTrans(:,2)) < 0.01);
%             toc
            
            %closest point to lidar accepted as return
%             tic
%             validIdxsNum = find(validIdxs);
%             toc
%             tic
            [range, returnIdx] = min(pointCloudTrans(validIdxs,3));
            
            if ~isempty(range)
                LidarImages(jj,kk,ii) = range;
            end
%             toc
        end
        disp(jj)
    end
end
toc

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



