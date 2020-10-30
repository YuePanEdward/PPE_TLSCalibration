%% Project Parameter Estimation (PPE) 2020 HS
% Main script
% Group: Yue & Josianne
%% Brief introduction
% Robust parameter estimation based on Gauss-Markov Model and Danish Method
% for terristial laser scanner intrinsic and extrinsic calibration
%% Denotations 
% unknown_vector x: [ APs, EP1, EP2, ... ]
% APs: [a0, b1, b2, c0],  the additional parameters for range(rou), horizontal angle(theta) 
% and elevation angle(alpha) correction 
% EPi: [omega, phi, kappa, t_x, t_y, t_z]_i, the extrinsic (geo-reference)
% parameters for each scan i, which is defined according to the following
% equation:
% p_i   =  R_io * (p_o - t_oi)
% namely, 
% p_i   =  Rx(omega) * Ry(phi) * Rz(keppa) * (p_o - [t_x, t_y, t_z]'),
% where p_i is the point coordinate in scanner i's coordinate system
% p_o is the point coordinate in external object coordinate system


%% Preparation
close all
clear all
clc

disp('PPE TLS calibration solution');
disp('By Yue Pan & Fandré Josianne');
% Define constants
% measurment standard deviation:
deg2rad_ratio=pi/180;
sigma_0 = 0.001;                       % prior unit weight standard deviation (m)
sigma_rho = 0.002;                     % range measurment (m)
sigma_theta = 0.005*deg2rad_ratio;     % horizontal angle (rad) 
sigma_alpha = 0.005*deg2rad_ratio;     % vertical angle (rad)

% delta_d for numerical derivative
delta_t = 1e-8;                        

% iteration thresholds
max_in_iter=20;                        % max iteration number for the internal loop (for gauss-markov parameter estimation)
max_ex_iter=10;                        % max iteration number for the external loop (for danish outlier detection and removal)
incre_ratio_thre = 1e-5;               % internal loop convergence condition
danish_converge_thre = 1e-2;           % external loop convergence condition

% prior knowledge for observation plausibility check
max_rou = 10.0;                        % max pausible range measurement value (m)  
max_alpha = 80*deg2rad_ratio;          % max pausible elevation angle value (rad)

% outlier threshold for danish method
danish_ratio_thre = 2.0;               % x sigma

%% I. Import data

disp('---------I. Data Import--------');

% Set the data path and prefix here
% [Test_data]
% data_path = strcat ('Test_data', filesep, 'Testdata_1');        
% data_prefix = 'PPE_TLS_T1_';  
% data_path = strcat ('Test_data', filesep, 'Testdata_2');      
% data_prefix = 'PPE_TLS_T2_';   

% [Final_data (with outlier)]
data_path = strcat ('Final_data', filesep, 'Finaldata_1');        
data_prefix = 'PPE_TLS_Fa_';  
% data_path = strcat ('Final_data', filesep, 'Finaldata_2');      
% data_prefix = 'PPE_TLS_Fb_'; 

% use filesep to represent / or \ on different operating system
scans_path=strcat(data_path, filesep, 'Scans', filesep); 
ops_path=strcat(data_path, filesep, 'OPs', filesep, data_prefix, 'OP.txt'); 
% import all the object points(OPs) in object's cartesian coordinate system
ops_raw=Read_OPFile(ops_path); %import the struct for ops
ops=ops_raw.XYZ';  % N * 3 matrix
op_count=size(ops,1); % N

disp(['Import [',num2str(op_count) ,'] OPs from [', ops_path , '] done']);

% import all the scan measurements in scanner's cartesian coordinte system
scans_in_cart = Read_ScanFile(scans_path, data_prefix); 
scan_count=length(scans_in_cart); % s
disp(['Import [', num2str(scan_count), '] Scans from [', scans_path ,'] done']);

%% II. Convert the scan measurements from Cartesian to spherical coordinate system

disp('---------II. Cart2Sphe--------');

% allocate
scans_in_sphe=cell(1, scan_count);

for i=1:scan_count
    scans_in_sphe{i} = cart2sphe(scans_in_cart{i}); % x,y,z --> rou,alpha,theta (unit:rad)
end

disp(['Cart2Spher convertion done for [',num2str(scan_count), '] scans']);

%% (Optional: mannually add some outliers)



%% III. (Optional) plausibility check
disp('---------III. Pre-plausibility check--------');

point_outlier_mask=ones(op_count * scan_count,1); % take each point
pre_outlier_count=0;

% Regard those measurements with too large range and elevation angle as outlier
for i=1:scan_count  
    for j=1:op_count
        temp_rou = scans_in_sphe{i}(j,1);
        temp_alpha = scans_in_sphe{i}(j,3);
        if(temp_rou > max_rou || temp_alpha > max_alpha)    
            point_outlier_mask((i-1)*op_count+j)=0; % the point (3 observations) would be regard as outlier
            pre_outlier_count=pre_outlier_count+1;  
        end 
    end
end
disp(['Pre-plausibility check done for [',num2str(op_count * scan_count), '] measurements, [', num2str(pre_outlier_count), '] outliers found.']);   

outlier_mask_rep = repmat(point_outlier_mask, [1,3])';
outlier_mask= outlier_mask_rep(:); % take each coordinate measurement

%% IV. Get initial guess

disp('---------IV. Assign initial guess of unknowns--------');

% initial guess for unknown vector (APs, EPs): (4+6s x 1)
ap_count= 4;
unknown_count = ap_count+6*scan_count; % u
% allocate
x_0=zeros(unknown_count, 1); % the initial value for 4 APs are assigned as 0

% for each scan, estimate the initial guess of eps (To,s: the 6DOF transformation from scanner to objective coordinate system)
for i=1:scan_count
   % input: source , target point cloud (matched with the same order),
   % output: transformation parameters: [omega(rx), phi(ry), kappa(rz), tx, ty, tz] %
   % tran source (station) to target (object)
   
   % only use the measurements pass the pre-plausibility check
   % (point_outlier_mask)
   x_0(ap_count+6*(i-1)+1:ap_count+6*(i-1)+6)=EstimateTranFromCorr(scans_in_cart{i}, ops); % using SVD (Horn, 1987)
end

disp(['Assign the initial value for [', num2str(unknown_count), '] unknowns done']);
disp_unknown_vector(x_0,ap_count,scan_count);

%% V. Conduct adjustment

disp('---------V. Apply adjustment iteratively--------');

ob_count = 3 * op_count * scan_count; % n = 3Ns

% assign observation vector y
y = zeros(ob_count,1); 
for i=1:scan_count
   y((i-1)*3*op_count+1:i*3*op_count)= reshape(scans_in_sphe{i}',3*op_count,1); % in spherical coordinate
end

disp(['Assign the obsevartion vector for [', num2str(ob_count), '] observations done']);

% stochastic model (K, sigma_0, Q, P)
% construct variance-covariance matrix
sigma_obs=[sigma_rho,sigma_theta,sigma_alpha];
sigma_obs_repeat=repmat(sigma_obs',ob_count/3,1);
D_mat=diag(sigma_obs_repeat.^2);  % variance-covariance matrix D
% sigma_0: a prior unit weight standard deviation              
Q_mat = (1.0/(sigma_0^2))*D_mat;  % cofactor matrix
P_mat_0 = inv(Q_mat);             % initial weight matrix

disp('Assign the initial weight matrix done');

x_temp=x_0; % assign initial value for unknown vector
P_mat_temp=P_mat_0;
ex_iter_count=0;

% Introduce the structure temp_adjustment_data
% the definition of each field: 
%  - x: initial guess of unknown vector
%  - y: observation vector
%  - P: weight matrix
%  - sigma_0: unit weight standard deviation
%  - dt: auto-derivative increment value
%  - op: object points' coordinate in cartesian system
%  - scans: a cell storing the measurements of the object points in each
%  scan's spherical system
%  - ap_count: number of the additional parameters
%  - max_iter_count: maximum iteration number of the inner-loop
%  - incre_ratio_thre: threshold for the ratio of the unknown value's
%  increment of the adjacent iterations, once the largest increment ratio
%  among the elements in the unknown vector is smaller than this threshold,
%  the loop would be regarded as converged
%  - outlier_mask: mask vector for the measurements, 0 and 1 indicate the
%  outlier and inlier respectively
temp_adjustment_data = struct('x',x_temp, 'y',y, 'P',P_mat_temp, 'sigma_0',...
   sigma_0, 'dt', delta_t, 'op', ops, 'scans', {scans_in_sphe}, 'ap_count',...
   ap_count, 'max_iter_count', max_in_iter, 'incre_ratio_thre', incre_ratio_thre,...
   'outlier_mask', outlier_mask);
 
% [solved issue]: to put the cell array into a struct, you need to use { }
% outside the cell array to make it a cell

iter_count=1; % external loop iteration number
danish_converged=0;
w_vec_last=diag(P_mat_0);

% External loop for detecting and removing outliers (Danish method)
while (~danish_converged && iter_count < max_ex_iter)
% iterate until reaching the termination criteria
   
   disp(['--------Danish Method iteration [', num2str(iter_count), ']--------']);
   % Adjustment calculation
   [x_p, Q_xx_mat, res_vec]= RunGMMAdjust(temp_adjustment_data);
   
   % Danish method weight updating
   % TODO: mask those measurements with too small weight as outliers
   % find the outlier index among current inliers according to dannish method
   danish_index = find(temp_adjustment_data.outlier_mask.*abs(res_vec)./sigma_obs_repeat > danish_ratio_thre); 
   danish_ratio = ones(ob_count,1);
   
   % Danish method coefficient
   if (iter_count<3) 
          danish_f=4.4;
   else
          danish_f=3.0;
   end
   danish_ratio(danish_index) = (exp(-(abs(res_vec(danish_index))./sigma_obs_repeat(danish_index)).^danish_f)).^0.05;
   w_vec = diag(P_mat_0).* danish_ratio; % diagonal weight vector
   temp_adjustment_data.P = diag(w_vec); % update weight matrix
   
   disp(['[' , num2str(size(danish_index,1)), '] candidate outliers found by Danish Method, downweight them.']);
   
   % Judge convergence
   d_w = abs(w_vec-w_vec_last)./diag(P_mat_0);
   if(max(d_w) < danish_converge_thre)
       danish_converged=1;
       disp('Danish method loop converged');
   end 

   iter_count=iter_count+1;
   w_vec_last=w_vec;
   
end

%% VI. Results output and visualization

disp('---------VI. Calibration results output and visualization--------');

% statistics
D_xx=sigma_0^2 * Q_xx_mat;                   % variance-covariance matrix of the estimated unknowns
sigma_xx=sqrt(diag(D_xx));                   % standard deviation of the estimated unknowns 
C_xx=D_xx./(sigma_xx*sigma_xx');             % correlation matrix of the estimated unknowns

% final unkonown estimation
disp('Final adjustment results of the unknown vector:');
disp_unknown_vector_std(x_p,sigma_xx,ap_count,scan_count);

% final outliers
outliers_by_plausibility_check=find(~outlier_mask);
outliers_by_danish_method=danish_index;
outliers_by_plausibility_check_count = size(outliers_by_plausibility_check,1);
outliers_by_danish_method_count = size(outliers_by_danish_method,1);
total_outlier_count=outliers_by_plausibility_check_count+outliers_by_danish_method_count;
fprintf('[%d] outliers found, [%d] by plausibility check, [%d] by danish method\n', total_outlier_count,outliers_by_plausibility_check_count,outliers_by_danish_method_count);

observation_status=zeros(ob_count,1);
observation_status(outliers_by_plausibility_check)=1;
observation_status(outliers_by_danish_method)=2;
measurement_status=cell(1,scan_count);
for i=1:scan_count
    scan_obs=observation_status((i-1)*3*op_count+1: i*3*op_count);
    measurement_status{i}=reshape(scan_obs,[3,op_count])';
end


unknown_name={'a0','b1','b2','c0'};
for i=1:scan_count
   unknown_name = {unknown_name{:}, ['omega', num2str(i)],['phi', num2str(i)],['kappa', num2str(i)],['X', num2str(i)],['Y', num2str(i)],['Z', num2str(i)]};
end

% Plot the posterior covariance matrix
figure(1);
set(gcf,'Position',[100 100 1800 600])
subplot(121);
heatmap(unknown_name, unknown_name, abs(D_xx), 'ColorScaling','log', 'Colormap', jet, 'Title', 'Posterior Covariance');

% Plot the posterior correlation matrix
subplot(122);
heatmap(unknown_name, unknown_name, C_xx, 'Colormap', jet, 'Title', 'Posterior Correlation');

% Plot the scanner and the ops
figure(2);
for i=1:scan_count
   plot_scanner(x_p(i*6-1:i*6+4),i,2);
end
scatter3(ops(:,1),ops(:,2),ops(:,3),15,'m','filled');
xlabel('X(m)');
ylabel('Y(m)');
zlabel('Z(m)');
%legend('scanner position', 'scanner x-axis','scanner y-axis', 'scanner z-axis', 'ops', 'Location','southoutside');
title('Overview of the scanners and the points');

% Plot the measurement of each scanner
figure(3);
for i=1:scan_count
   subplot(1,scan_count,i);
   plot_scanner(zeros(6,1),i,3);
   inliers=scans_in_cart{i}((max(measurement_status{i}'))==0,:);
   pausibility_outliers=scans_in_cart{i}((max(measurement_status{i}'))==1,:);
   dannish_outliers=scans_in_cart{i}((max(measurement_status{i}'))==2,:);
   scatter3(inliers(:,1),inliers(:,2),inliers(:,3),15,'g','filled');
   scatter3(pausibility_outliers(:,1),pausibility_outliers(:,2),pausibility_outliers(:,3),15,'r','filled');
   scatter3(dannish_outliers(:,1),dannish_outliers(:,2),dannish_outliers(:,3),15,'b','filled');
   xlabel('X(m)');
   ylabel('Y(m)');
   zlabel('Z(m)');
   legend('scanner position', 'scanner x-axis','scanner y-axis', 'scanner z-axis', 'inlier measurement', 'outlier detected by plausibility check', 'outlier detected by Danish method','Location','southoutside');
end
%title('measurements of each scanner (green: inlier, red: outlier from plausibility check, blue: outlier from danish method)');




