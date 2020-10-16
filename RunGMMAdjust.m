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

% Pre-Allocation
iter_count = 0;
d_xi = inf;
crit = 1;

% Iteration % TODO:
while crit && (iter_count < adjustment_data_struct.max_iter_count)
    iter_count = iter_count + 1;
    
    % Caculate A matrix (Jocubians) by numerical derivative 
    % Allocate
    A_mat = zeros(ob_count,unknown_count);
    for i = 1:scan_count
       for j = 1:op_count
          for k = 1:ap_count % calculate derivative for APs
              ObsFunction(op, sp, aps, eps)
               
               A(k:k+2,k) = derivative(@funcModel,q,xi(1:4),xi(j*6-1:j*6+4),OP{j}(i,:),scans{j}(i,:),dt(q));
          end

           for q = (4+ 6*j-5):(4+ 6*j) % calculate derivative for EPs
               A(k:k+2,q) = derivative(@funcModel,(q-6*(j-1)),xi(1:4),xi(j*6-1:j*6+4),OP{j}(i,:),scans{j}(i,:),dt(q));
           end
           k = k + 3;
    end
    
    
    A = zeros(3*sum([n{:}]),4+6*s);
    k = 1;
    for j = 1:s
       for i = 1:n{j}
           for q = 1:4
               A(k:k+2,q) = derivative(@funcModel,q,xi(1:4),xi(j*6-1:j*6+4),OP{j}(i,:),scans{j}(i,:),dt(q));
           end

           for q = (4+ 6*j-5):(4+ 6*j)
               A(k:k+2,q) = derivative(@funcModel,(q-6*(j-1)),xi(1:4),xi(j*6-1:j*6+4),OP{j}(i,:),scans{j}(i,:),dt(q));
           end
           k = k + 3;
       end

    end

end
