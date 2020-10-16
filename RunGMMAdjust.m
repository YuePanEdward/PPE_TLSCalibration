function [x_p, Q_xx_mat, res_vec] = RunGMMAdjust(adjustment_data_struct)
%RunGMMAdjust: do iterative incremental parameter estimation based on the Gauss-Markow Model(GMM)         

x_0=adjustment_data_struct.x;
y=adjustment_data_struct.y;
P=adjustment_data_struct.P;
sigma_0=adjustment_data_struct.sigma_0;
dt=adjustment_data_struct.dt;
ops=adjustment_data_struct.op;
scans_sphe=adjustment_data_struct.scans;
ap_count=adjustment_data_struct.ap_count;

scan_count=length(scans_sphe);  %s
op_count=size(ops,1);           %N
unknown_count=ap_count+6*scan_count;   %u
ob_count=3*op_count*can_count;  %n


iter_count = 0;
is_converged =0;

% Iteration % TODO:
while (iter_count < adjustment_data_struct.max_iter_count)
    iter_count = iter_count + 1;
    
    % Caculate A matrix (Jocubians) by numerical derivative 
    % Allocate
    A_mat = zeros(ob_count,unknown_count);
    for i = 1:scan_count
       for j = 1:op_count
          for k = 1:ap_count % calculate derivative for APs
              cur_res = ObsFunction(op, sp, aps, eps);
              A_mat() =  Derivative(cur_res, dt);
      
          end
          for k = (6*i-1):(6*i+4) % calculate derivative for EPs
              cur_res = ObsFunction(op, sp, aps, eps);
              A_mat() =  Derivative(cur_res, dt);
          end
       end
    end
    %TODO
    

end
