% Generates Simulated Lidar measurements from the viewpoint of cam0 in the
% frame of cam0
% Corey Marcus
% UT Austin: Aerospace Engineering

clear
close all
clc

%% Options
LidarFOVHeight = pi/6; %radians
LidarFOVWidth = pi/4; %radians
% LidarFOVHeight = .01; %radians
% LidarFOVWidth = .01; %radians
LidarArrayWidth = 5;
LidarArrayHeight = 5;
lidarPointSwath = pi/180; %1 degree lidar point return
lidarImageRate = 1; %one to create lidar image for every camera image
trackPoints = true; %tracks points pinged for debugging

%% Main

%maximum lidar FOV, used to restrict point cloud searching
lidarFOVMax = max([LidarFOVHeight/2, LidarFOVWidth/2]);

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
LidarAzAngles = linspace(-LidarFOVWidth/2, LidarFOVWidth/2,LidarArrayWidth); 
LidarElAngles = linspace(LidarFOVHeight/2, -LidarFOVHeight/2,LidarArrayHeight); %note intentional sign reversal

%create a calibration output
FID = fopen('lidarCalib.csv','w');
fprintf(FID,'%3i, %3i, \n', LidarArrayWidth, LidarArrayHeight);
for ii = 1:LidarArrayWidth
    fprintf(FID,'%8f, ',LidarAzAngles(ii));
end
fprintf(FID,'\n');
for ii = 1:LidarArrayHeight
    fprintf(FID,'%8f, ',LidarElAngles(ii));
end
fclose(FID);

%create a lidar image array
LidarImages = zeros(LidarArrayHeight,LidarArrayWidth,L);

%initialize return point array for debugging
returnPtMat = [];

%Create a quaternion to rotate between a camera frame and an azimuth,
%elevation frame
Q_lidar2azEl = angle2quat(-pi/2,pi/2,0,'YXZ');

%images to operate on
% imageIdx = 1:lidarImageRate:L;
imageIdx = 1;

parfor ii = imageIdx %iteration on images
    
    %shift point cloud to camera origin
    pointCloudShift = pointCloud - P_inertial2cam(:,matchedInertialIdxs(ii))';
    
    %convert point cloud to spherical coordinates centered on camera
    [az, el, r] = cart2sph(pointCloudShift(:,1),pointCloudShift(:,2),pointCloudShift(:,3));
    
    %create euler angles for current lidar bore sight
    [boreYaw, borePitch, boreRoll] = quat2angle(...
        quatmultiply(Q_inertial2cam(:,matchedInertialIdxs(ii))',...
        Q_lidar2azEl),'ZYX');
    
    %yaw = azimuth, but pitch = negative elevation
    boreAz = boreYaw;
    boreEl = -borePitch;
    
    %purge points far out of view of lidar
    boreAzMin = boreAz - lidarFOVMax - pi/6;
    boreAzMax = boreAz + lidarFOVMax + pi/6;
    boreElMin = boreEl - lidarFOVMax - pi/6;
    boreElMax = boreEl + lidarFOVMax + pi/6;
    
    %different cases depending on if yaw and pitch need wrapping
    if (boreAzMin > -pi) && (boreAzMax < pi) %nominal case
        azValid = (az > boreAzMin) & (az < boreAzMax);
    elseif (boreAzMin <= -pi) && (boreAzMax < pi) %wrap below
        boreAzMin = boreAzMin + 2*pi;
        azValid = (az > boreAzMin) | (az < boreAzMax);
    elseif (boreAzMin > -pi) && (boreAzMax >= pi) %wrap above
        boreAzMax = boreAzMax - 2*pi;
        azValid = (az > boreAzMin) | (az < boreAzMax);
    else
        disp('Bounding Error!')
    end
    
    if (boreElMin > -pi/2) && (boreElMax < pi/2) %nominal case
        elValid = (el > boreElMin) & (el < boreElMax);
    elseif (boreElMin <= -pi/2) && (boreElMax < pi/2) %wrap below
        boreElMin = boreElMin + pi;
        elValid = (el > boreElMin) | (el < boreElMax);
    elseif (boreElMin > -pi/2) && (boreElMax >= pi/2) %wrap above
        boreElMax = boreElMax - pi;
        elValid = (el > boreElMin) | (el < boreElMax);
    elseif (boreElMin <= -pi/2) && (boreElMax >= pi/2) %all points visible
        elValid = true(length(el),1);
    else
        disp('Bounding Error!')
    end
    
    %indicies within range of lidar
    validIdx = azValid & elValid;
    
    %valid points
    azValid = az(validIdx);
    elValid = el(validIdx);
    rvalid = r(validIdx);
    
    for jj = 1:LidarArrayHeight
        for kk = 1:LidarArrayWidth
            
            %create quaternion corresponding to this angle
            Q_cam2lidarPt = angle2quatComp(LidarAzAngles(kk),LidarElAngles(jj),0,'ZYX');
            
            %create quaternion from inertial to this array point
            Q_inertial2lidarPt = quatmultiply(...
                quatmultiply(Q_inertial2cam(:,matchedInertialIdxs(ii))',...
                Q_lidar2azEl), Q_cam2lidarPt);
            
            %create euler angles for current lidar point
            [ptYaw, ptPitch, ptRoll] = quat2angle(Q_inertial2lidarPt,'ZYX');
            
            %yaw = azimuth, but pitch = negative elevation
            ptEl = -ptPitch;
            ptAz = ptYaw;
            
            %bounds for lidar point swath
            ptAzMin = ptAz - lidarPointSwath/2;
            ptAzMax = ptAz + lidarPointSwath/2;
            ptElMin = ptEl - lidarPointSwath/2;
            ptElMax = ptEl + lidarPointSwath/2;
            
            %find vaild points
            %different cases depending on if yaw and pitch need wrapping
            if (ptAzMin > -pi) && (ptAzMax < pi) %nominal case
                azValidPt = (azValid > ptAzMin) & (azValid < ptAzMax);
            elseif (ptAzMin <= -pi) && (ptAzMax < pi) %wrap below
                ptAzMin = ptAzMin + 2*pi;
                azValidPt = (azValid > ptAzMin) | (azValid < ptAzMax);
            elseif (ptAzMin > -pi) && (ptAzMax >= pi) %wrap above
                ptAzMax = ptAzMax - 2*pi;
                azValidPt = (azValid > ptAzMin) | (azValid < ptAzMax);
            else
                disp('Bounding Error!')
            end
            
            if (ptElMin > -pi/2) && (ptElMax < pi/2) %nominal case
                elValidPt = (elValid > ptElMin) & (elValid < ptElMax);
            elseif (ptElMin <= -pi/2) && (ptElMax < pi/2) %wrap below
                ptElMin = ptElMin + pi;
                elValidPt = (elValid > ptElMin) | (elValid < ptElMax);
            elseif (ptElMin > -pi/2) && (ptElMax >= pi/2) %wrap above
                ptElMax = ptElMax - pi;
                elValidPt = (elValid > ptElMin) | (elValid < ptElMax);
            elseif (ptElMin <= -pi/2) && (ptElMax >= pi/2) %all points visible
                elValidPt = true(length(elValid),1);
            else
                disp('Bounding Error!')
            end
            
            %points within view of this particular lidar return
            validIdxPt = azValidPt & elValidPt;
            
            %extract candidate points
            candAz = azValid(validIdxPt);
            candEl = elValid(validIdxPt);
            candr = rvalid(validIdxPt);
            
            %assign range
            if ~isempty(candr)
                [ptRange, ptIdx] = min(candr);
                LidarImages(jj,kk,ii) = ptRange;
            end
            
            %track points for debugging
            if ~isempty(candr) && trackPoints
                [returnPtX, returnPtY, returnPtZ] = sph2cart(candAz(ptIdx),candEl(ptIdx),candr(ptIdx));
                returnPtMat = [returnPtMat;
                    [returnPtX, returnPtY, returnPtZ] + P_inertial2cam(:,matchedInertialIdxs(ii))'];
            end
            
        end
        
    end
    
    disp(ii/L);
    
end

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

% figure
% plot(P_inertial2body(1,:),P_inertial2body(2,:))
% hold on
% plot(P_inertial2cam(1,:),P_inertial2cam(2,:))
% scatter(P_inertial2body(1,1),P_inertial2body(2,1))
% xlabel('x')
% ylabel('y')
% quiver(P_inertial2cam(1,1:100:end),P_inertial2cam(2,1:100:end),...
%     Barrel_camInInertial(1,1:100:end),Barrel_camInInertial(2,1:100:end))
% axis equal
%
%
% figure
% plot3(P_inertial2body(1,:),P_inertial2body(2,:),P_inertial2body(3,:))
% hold on
% plot3(P_inertial2cam(1,:),P_inertial2cam(2,:),P_inertial2cam(3,:))
% scatter3(P_inertial2body(1,1),P_inertial2body(2,1),P_inertial2body(3,1))
% xlabel('x')
% ylabel('y')
% zlabel('z')
% quiver3(P_inertial2cam(1,1:100:end),P_inertial2cam(2,1:100:end),P_inertial2cam(3,1:100:end),...
%     Barrel_camInInertial(1,1:100:end),Barrel_camInInertial(2,1:100:end),Barrel_camInInertial(3,1:100:end))
% axis equal
%
% figure
% plot3(P_inertial2body(1,:),P_inertial2body(2,:),P_inertial2body(3,:))
% hold on
% plot3(P_inertial2cam(1,:),P_inertial2cam(2,:),P_inertial2cam(3,:))
% scatter3(P_inertial2body(1,1),P_inertial2body(2,1),P_inertial2body(3,1))
% scatter3(pointCloud(1:10:end,1),pointCloud(1:10:end,2),pointCloud(1:10:end,3))
% xlabel('x')
% ylabel('y')
% zlabel('z')
% quiver3(P_inertial2cam(1,1:100:end),P_inertial2cam(2,1:100:end),P_inertial2cam(3,1:100:end),...
%     Barrel_camInInertial(1,1:100:end),Barrel_camInInertial(2,1:100:end),Barrel_camInInertial(3,1:100:end),0)
% axis([-5 5 -5 5 0 3])

%plot camera position and return point positions
if trackPoints
    figure
    plot3(P_inertial2cam(1,:),P_inertial2cam(2,:),P_inertial2cam(3,:))
    hold on
    scatter3(returnPtMat(:,1), returnPtMat(:,2), returnPtMat(:,3))
    scatter3(P_inertial2body(1,1),P_inertial2body(2,1),P_inertial2body(3,1))
    quiver3(P_inertial2cam(1,matchedInertialIdxs(imageIdx)),...
        P_inertial2cam(2,matchedInertialIdxs(imageIdx)),...
        P_inertial2cam(3,matchedInertialIdxs(imageIdx)),...
        Barrel_camInInertial(1,matchedInertialIdxs(imageIdx)),...
        Barrel_camInInertial(2,matchedInertialIdxs(imageIdx)),...
        Barrel_camInInertial(3,matchedInertialIdxs(imageIdx)),0)
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    %plot point cloud and return point positions
    figure
    scatter3(pointCloud(1:1:end,1),pointCloud(1:1:end,2),pointCloud(1:1:end,3),...
    1,'filled')
    hold on
    scatter3(returnPtMat(:,1), returnPtMat(:,2), returnPtMat(:,3))
    scatter3(P_inertial2body(1,1),P_inertial2body(2,1),P_inertial2body(3,1))
    quiver3(P_inertial2cam(1,matchedInertialIdxs(imageIdx)),...
        P_inertial2cam(2,matchedInertialIdxs(imageIdx)),...
        P_inertial2cam(3,matchedInertialIdxs(imageIdx)),...
        Barrel_camInInertial(1,matchedInertialIdxs(imageIdx)),...
        Barrel_camInInertial(2,matchedInertialIdxs(imageIdx)),...
        Barrel_camInInertial(3,matchedInertialIdxs(imageIdx)),0)
    xlabel('x')
    ylabel('y')
    zlabel('z')
%     axis([-6 6 -5 5 -1 3])
    axis([-4.75 -4 .6 1.8 0.75 1])
end

