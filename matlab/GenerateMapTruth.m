% Generates truth map depth from the viewpoint of cam0 in the
% frame of cam0
% Corey Marcus
% UT Austin: Aerospace Engineering

clear
close all
clc

%% Options

% Specify the frames to generate truth for
% targs = 73:173;
targs = 73 + [0 35 54]; % These frames are often chosen as keyframes

% search swath
searchswath = 2*pi/180;

%% Main

% load dataset
datasetPath = '/home/cm58349/Documents/EuRoC_data/V2_01_easy'; % atlantis config
% datasetPath = 'C:\Users\cm58349\Documents\EuRoC_data\V2_01_easy';
% datasetPath =
% 'C:\Users\corey\Documents\SharedFolder\EuRoC_data\V2_01_easy'; % personal
% config
addpath('quaternion');
dataset = dataset_load(datasetPath);

% Camera calibration matrix
% This is the undistorted camera calibration model as reported by
% LSD-SLAM's openCV routines
K = [303.072       0 308.727
      0 418.033  250.18
      0       0       1];
Kinv = inv(K);

% Size of the cropped images used by LSD-SLAM
width = 640;
height = 480;

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

% Create truth images
Ntruth = length(targs);
truthimages = zeros(height,width,Ntruth);

% Begin looping
targidx = 1; % target index for parfor compatibility
for ii = targs
    %shift point cloud to camera origin
    pointCloudShift = pointCloud - P_inertial2cam(:,matchedInertialIdxs(ii))';
    
    %convert point cloud to spherical coordinates centered on camera
    [az, el, r] = cart2sph(pointCloudShift(:,1),pointCloudShift(:,2),pointCloudShift(:,3));

    % Parallelize with rows
    parfor jj = 1:height
        disp(jj)
        for kk = 1:width
            % Create a vector corresponding to this camera pixel
            v_cam = Kinv*[kk; jj; 1];
            v_cam = v_cam/norm(v_cam);

            % rotate vector into the inertial frame
            v_inert = quatrotate(quatconj(Q_inertial2cam(:,matchedInertialIdxs(ii))'),v_cam');

            % get euler angles for point
            ptelev = asin(v_inert(3));
            v_xy = sqrt(v_inert(1)^2 + v_inert(2)^2);
            ptazim = atan2(v_inert(2),v_inert(1));

            % bounds for lidar point swath
            ptazimmin = ptazim - searchswath/2;
            ptazimmax = ptazim + searchswath/2;
            ptelevmin = ptelev - searchswath/2;
            ptelevmax = ptelev + searchswath/2;

            % find vaild points
            % different cases depending on if yaw and pitch need wrapping
            if (ptazimmin > -pi) && (ptazimmax < pi) %nominal case
                azValidPt = (az > ptazimmin) & (az < ptazimmax);
            elseif (ptazimmin <= -pi) && (ptazimmax < pi) %wrap below
                ptazimmin = ptazimmin + 2*pi;
                azValidPt = (az > ptazimmin) | (az < ptazimmax);
            elseif (ptazimmin > -pi) && (ptazimmax >= pi) %wrap above
                ptazimmax = ptazimmax - 2*pi;
                azValidPt = (az > ptazimmin) | (az < ptazimmax);
            else
                warning('Bounding Error!')
            end
            
            if (ptelevmin > -pi/2) && (ptelevmax < pi/2) %nominal case
                elValidPt = (el > ptelevmin) & (el < ptelevmax);
            elseif (ptelevmin <= -pi/2) && (ptelevmax < pi/2) %wrap below
                ptelevmin = ptelevmin + pi;
                elValidPt = (el > ptelevmin) | (el < ptelevmax);
            elseif (ptelevmin > -pi/2) && (ptelevmax >= pi/2) %wrap above
                ptelevmax = ptelevmax - pi;
                elValidPt = (el > ptelevmin) | (el < ptelevmax);
            elseif (ptelevmin <= -pi/2) && (ptelevmax >= pi/2) %all points visible
                elValidPt = true(length(elValid),1);
            else
                warning('Bounding Error!')
            end
            
            % points within view of this particular lidar return
            valididxpt = azValidPt & elValidPt;
            
            % extract candidate points
            candr = r(valididxpt);

             %assign range
            if ~isempty(candr)
                
                %get point
                [ptrange, ~] = min(candr);
                
                %make sure it is not less than zero
                if(ptrange < 0)
                    ptrange = 0;
                    warning('Negative Range Found!')
                end
                
                %assign
                truthimages(jj,kk,targidx) = ptrange;
            else
                warning("no match")
            end
   
        end
    end
    
    % increment target index
    targidx = targidx + 1;
end

% Save data
save('MapTruthData.mat','truthimages','-v7.3');

%% Plotting
figure
mesh(truthimages(:,:,1))
