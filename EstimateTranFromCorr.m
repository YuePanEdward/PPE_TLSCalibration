function [tran_params] = EstimateTranFromCorr(source,target)
% EstimateTranFromCorr Function estimate the 6DOF rigid body transformation
% from two point cloud (source to target) with pre-determined correspondence
% Using singular value decomposition (SVD) closed-form solution
% refer to [Horn, 1987]

source_grav_cent=mean(source);
target_grav_cent=mean(target);
source_degrav=source-source_grav_cent;   % source (p)
target_degrav=target-target_grav_cent;   % target (q)
[U,S,V]=svd(source_degrav'*target_degrav); % apply SVD
R_mat=U*V'; % Rotation matrix
t_vec=target_grav_cent-(R_mat*source_grav_cent')'; % translation vector

% nav toolbox required for the rotm2eul function
r_vec=-rotm2eul(R_mat,'XYZ'); % roll, pitch, yaw (in rad , anti-clockwise)

tran_params=[r_vec t_vec];

end

