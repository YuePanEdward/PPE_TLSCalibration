function disp_unknown_vector_std(x_vec, sigma_x, ap_count, scan_count)
%disp_unknown_vector_std: display each element of the unknown vector and
%its standard deviation
% input: 
% - x_vec: unknown parameters (aps + eps of each scan)
% - sigma_x: unknown parameters' standard deviation
% - ap_count: number of aps
% - scan_count: number of scans

deg2rad_ratio=pi/180;

% disp('APs:');
% for i=1:ap_count
%     fprintf('%.3f\t',x_vec(i)); 
% end
% fprintf('\n'); 
disp('APs:')
fprintf('a0 = %8.3f ( mm )\t+-%5.2f ( mm )\n',x_vec(1)*1e3, sigma_x(1)*1e3);
fprintf('b1 = %8.3f (mdeg)\t+-%5.2f (mdeg)\n',x_vec(2)/deg2rad_ratio*1e3, sigma_x(2)/deg2rad_ratio*1e3);
fprintf('b2 = %8.3f (mdeg)\t+-%5.2f (mdeg)\n',x_vec(3)/deg2rad_ratio*1e3, sigma_x(3)/deg2rad_ratio*1e3);
fprintf('c0 = %8.3f (mdeg)\t+-%5.2f (mdeg)\n\n',x_vec(4)/deg2rad_ratio*1e3, sigma_x(4)/deg2rad_ratio*1e3);

for i=1:scan_count
    disp(['EPs of scan [', num2str(i), ']:']); 
    fprintf('omega = %8.3f (deg)\t+-%5.2f (mdeg)\n',x_vec(ap_count+6*i-5)/deg2rad_ratio, sigma_x(ap_count+6*i-5)/deg2rad_ratio*1e3);
    fprintf('phi   = %8.3f (deg)\t+-%5.2f (mdeg)\n',x_vec(ap_count+6*i-4)/deg2rad_ratio, sigma_x(ap_count+6*i-4)/deg2rad_ratio*1e3);
    fprintf('kappa = %8.3f (deg)\t+-%5.2f (mdeg)\n',x_vec(ap_count+6*i-3)/deg2rad_ratio, sigma_x(ap_count+6*i-3)/deg2rad_ratio*1e3);
    fprintf('X_s   = %8.3f ( m )\t+-%5.2f ( mm )\n',x_vec(ap_count+6*i-2), sigma_x(ap_count+6*i-2)*1e3);
    fprintf('Y_s   = %8.3f ( m )\t+-%5.2f ( mm )\n',x_vec(ap_count+6*i-1), sigma_x(ap_count+6*i-1)*1e3);
    fprintf('Z_s   = %8.3f ( m )\t+-%5.2f ( mm )\n\n',x_vec(ap_count+6*i), sigma_x(ap_count+6*i)*1e3);
end

end
