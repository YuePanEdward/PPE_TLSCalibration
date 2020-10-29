function disp_unknown_vector(x_vec, ap_count, scan_count)
%disp_unknown_vector: display each element of the unknown vector

disp('APs:');
for i=1:ap_count
    fprintf('%f\t',x_vec(i)); 
end
fprintf('\n'); 

for i=1:scan_count
    disp(['EP of scan [', num2str(i), ']:']); 
    fprintf('omega(rad): %f\t phi(rad): %f\t kappa(rad): %f\n',x_vec(ap_count+6*i-5:ap_count+6*i-3)); 
    fprintf('X(m): %f\t Y(m): %f\t Z(m): %f\n',x_vec(ap_count+6*i-2:ap_count+6*i)); 
end

end

