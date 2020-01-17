function [LidarImages] = lidarImageBuild(LidarArrayHeight,...
    LidarArrayWidth, LidarYawAngles, LidarPitchAngles, Q_inertial2cam,...
    matchedInertialIdxs, P_inertial2cam, pointCloud)
%creates psuedo lidar range images

LidarImages = zeros(LidarArrayHeight,LidarArrayWidth,1);

for ii = 1:1 %iteration on images
    for jj = 1:LidarArrayHeight
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
            [range, ~] = min(pointCloudTrans(validIdxs,3));
            
            if ~isempty(range)
                LidarImages(jj,kk,ii) = range;
            end
%             toc
        end
        disp(jj)
    end
end

end

