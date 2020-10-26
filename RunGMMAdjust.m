function [x_p, Q_xx_mat, res_vec] = RunGMMAdjust(adjustment_data_struct)
%RunGMMAdjust: do iterative incremental parameter estimation based on the Gauss-Markow Model(GMM)         

x_0=adjustment_data_struct.x;
y=adjustment_data_struct.y;
P_mat=adjustment_data_struct.P;
sigma_0=adjustment_data_struct.sigma_0;
dt=adjustment_data_struct.dt;
ops=adjustment_data_struct.op;
scans_sphe=adjustment_data_struct.scans;
ap_count=adjustment_data_struct.ap_count;

scan_count=length(scans_sphe);  %s
op_count=size(ops,1);           %N
unknown_count=ap_count+6*scan_count;   %u
ob_count=3*op_count*scan_count;  %n


iter_count = 1;
is_converged =0;

% Internal Iteration 
% TODO
while (iter_count < adjustment_data_struct.max_iter_count)
     
    % Caculate A matrix (Jocubians) by numerical derivative 
    % Allocate
    ob_index = 1;
    A_mat = zeros(ob_count,unknown_count);
    b_vec = zeros(ob_count,1);
    
    for i = 1:scan_count
       cur_scan_sphe = scans_sphe{1,i};
       for j = 1:op_count
          for k = 1:unknown_count % calculate derivative for APs
              A_mat(ob_index:ob_index+2, k) = Derivative(ops(j,:), cur_scan_sphe(j,:), x_0(1:ap_count), x_0(ap_count+(i-1)*6+1:ap_count+i*6),dt, k);
              
              funct= ObsFunction(ops(j,:), cur_scan_sphe(j,:), x_0(1:ap_count), x_0(ap_count+(i-1)*6+1:ap_count+i*6));
              y_cur =y(ob_index:ob_index+2,1);
              b_vec(ob_index:ob_index+2 ,1) = y_cur - funct;
          end
          ob_index = ob_index+3;
       end
    end
    
    
    % Conduct adjustment
    x_i =inv( A_mat' * A_mat ) * (A_mat' * b_vec);
    
    %disp(['Increment of the unknown vector of iteration [', num2str(iter_count) , ']:', num2str(x_i)]);
    
    x_0 = x_0 + x_i;
    iter_count = iter_count + 1;
    
end

x_p = x_0; % why do you not give out just x_0 instead of saving it separatly as x_p
%disp(['Adjustment results of the unknown vector: ', num2str(x_i)]);

% TODO
Q_xx_mat=0;
res_vec=0; 


