function pol_points= cart2sphe(cart_points)
% cart2sphe: transformation from cartesian to polar (spherical) coordinates system.
%   Input:  xyz (n x 3) 'cart_points' in cartesian coordinate system
%   Output: rho theta alpha (n x 3) 'pol_points' in polar coordinate
%   system, angle unit: rad

rho = batchnorm(cart_points);                                  % range
theta = atan2(cart_points(:,2),cart_points(:,1));               % horizontal angle (rad)
alpha = atan2(cart_points(:,3),batchnorm(cart_points(:,1:2))); % vertical (elevation) angle (rad)
pol_points= [rho,theta,alpha];

end