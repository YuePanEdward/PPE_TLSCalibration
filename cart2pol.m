function pol_points= cart2pol(cart_points)
% cart2pol: transformation from cartesian to polar (spherical) coordinates system.
%   Input:  xyz (n x 3) 'cart_points' in cartesian coordinate system
%   Output: rho theta alpha (n x 3) 'pol_points' in polar coordinate system 

rho = vectornorm(cart_points);                                  % range
theta = atan2(cart_points(:,2),cart_points(:,1));               % horizontal angle
alpha = atan2(cart_points(:,3),vectornorm(cart_points(:,1:2))); % vertical (elevation) angle
pol_points= [rho,theta,alpha];

end

