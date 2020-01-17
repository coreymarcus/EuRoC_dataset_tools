function [LidarImages] = lidarImageBuild(LidarArrayHeight,...
    LidarArrayWidth, LidarYawAngles, LidarPitchAngles, NumImages,...
    Q_inertial2cam, matchedInertialIdxs, P_inertial2cam, pointCloud)
%creates psuedo lidar range images

LidarImages = zeros(LidarArrayHeight,LidarArrayWidth,NumImages);

for ii = 1:NumImages %iteration on images
    parfor jj = 1:LidarArrayHeight
        for kk = 1:LidarArrayWidth
            %             tic
            %create quaternion corresponding to this angle
            Q_cam2lidarray = angle2quatComp(LidarYawAngles(kk),LidarPitchAngles(jj),0,'ZXY');
            
            %generate a lidar ray in the inertial frame
            Q_inertial2lidarray = quatmultiply(Q_inertial2cam(:,matchedInertialIdxs(ii))',...
                Q_cam2lidarray);
            
            %transform point cloud to lidar frame
            shiftMat = repmat(P_inertial2cam(:,matchedInertialIdxs(ii))',size(pointCloud,1),1);
            pointCloudTrans = quatrotateComp(Q_inertial2lidarray, ...
                pointCloud - shiftMat);
            %             toc
            
            %find valid indicies in the point cloud trans
%             tic
            validPointCloud = pointCloudTrans((pointCloudTrans(:,3) > .1) & (abs(pointCloudTrans(:,1)) < 0.01) ...
                & (abs(pointCloudTrans(:,2)) < 0.01),:);
            
%             toc

            %closest point to lidar accepted as return
            %             tic
            %             validIdxsNum = find(validIdxs);
            %             toc
            %             tic
            
            
            if ~isempty(validPointCloud) > 0
                %                 [range, ~] = min(pointCloudTrans(validIdxs,3));
                [range, ~] = min(validPointCloud(:,3));
                LidarImages(jj,kk,ii) = range;
            end
            %             toc
        end
        %         disp(jj)
    end
    disp(ii/NumImages)
end

end

