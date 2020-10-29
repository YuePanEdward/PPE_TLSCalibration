function obs = ObsFunction(op, sp, aps, eps)
    %ObsFunction: the functional model for the GMM parameter estimation
    % Calculated in spherical coordinate (rho, theta, alpha)]
    
    % input: op (1 x 3)  X,Y,Z in cardestian coordinate system
    %        sp (1 x 3)  
 
    
    alpha = sp(3);
    d_rho = aps(1);
	d_theta = aps(2)*sec(alpha) + aps(3)*tan(alpha);
	d_alpha = aps(4);
    d_spher_coor = [d_rho, d_theta, d_alpha]; % 1 x 3
    
    % nav toolbox required for eul2rotm function
    % R_so
    R_mat=eul2rotm(eps(1:3)','XYZ'); % reconstruct rotation matrix from euler angles (anti-clockwise,rad)
    % t_os
    t_vec = eps(4:6); % translation vector
    
    % p_s = R_so * (p_o-t_os)
    est_scan_in_cart = (R_mat * (op'-t_vec))';  % 1 x 3
    
    est_scan_in_sphe = cart2sphe(est_scan_in_cart); % 1 x 3
	
    % get the function defined in spherical coordinate
	obs = (est_scan_in_sphe + d_spher_coor)'; % 3 x 1
    
end