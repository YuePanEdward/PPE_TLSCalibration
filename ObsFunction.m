function res = ObsFunction(op, sp, aps, eps)
    %ObsFunction: the functional model for the GMM parameter estimation
    % Calculated in spherical coordinate (rho, theta, alpha)
    
    alpha = sp(3);
    d_rho = aps(1);
	d_theta = aps(2)*sec(alpha) + aps(3)*tan(alpha);
	d_alpha = aps(4);
    d_spher_coor = [d_rho, d_theta, d_alpha];
    
    % nav toolbox required for eul2rotm function
    R_mat=eul2rotm(-eps(1:3)','XYZ'); % reconstruct rotation matrix from euler angles (anti-clockwise,rad)
    t_vec = eps(4:6); % translation vector
    
    scan_in_cart = (R_mat * (op'-t_vec))'; 
	
    % get resiudal defined in spherical coordinate
	res = d_spher_coor + cart2sphe(scan_in_cart);
    
end