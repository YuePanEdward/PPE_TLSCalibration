%% Project Parameter Estimation (PPE) 2020 HS
% Main script
% Group: Yue & Josianne
%% Brief introduction
% Robust parameter estimation based on Gauss-Markov Model and Danish Method
% for terristial laser scanner intrinsic and extrinsic calibration
%% Denotations


%% Preparation
close all
clear all
clc

disp('PPE TLS calibration software');

% Define constants
% measurment standard deviation:
%(guess)
deg2rad_ratio=pi/180;
sigma_0 = 0.001;                      % prior unit weight standard deviation (m)
sigma_rho = 0.005;                    % range measurment (m)
sigma_theta = 0.01*deg2rad_ratio;     % horizontal angle (rad) 
sigma_alpha = 0.01*deg2rad_ratio;     % vertical angle (rad)

% iteration thresholds
max_in_iter=20;
max_ex_iter=1;
incre_ratio_thre = 1e-3;


%% I. Import data

disp('---------I. Data Import--------');

% Set the data path and prefix here
data_path = 'Testdata_1';        
data_prefix = 'PPE_TLS_T1_';  
% data_path = 'Testdata_2';        
% data_prefix = 'PPE_TLS_T2_';   

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

%% III. Get initial guess

disp('---------III. Assign initial guess of unknowns--------');

% initial guess for unknown vector (APs, EPs): (4+6s x 1)
ap_count= 4;
unknown_count = ap_count+6*scan_count; % u
% allocate
x_0=zeros(unknown_count, 1); % the initial value for 4 APs are assigned as 0

% for each scan, estimate the initial guess of eps (To,s: the 6DOF transformation from scanner to objective coordinate system)
for i=1:scan_count
   % input: source , target point cloud (matched with the same order),
   % output: transformation parameters: [roll, pitch, yaw, tx, ty, tz] %
   % tran source (station) to target (object)
   x_0(ap_count+6*(i-1)+1:ap_count+6*(i-1)+6)=EstimateTranFromCorr(scans_in_cart{i}, ops); % using SVD (Horn, 1987)
end

disp(['Assign the initial value for [', num2str(unknown_count), '] unknowns done']);

%% IV. Conduct adjustment

disp('---------IV. Apply adjustment iteratively--------');

delta_t = 1e-8; % delta_d for numerical derivative

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
K_mat=diag(sigma_obs_repeat);  % variance-covariance matrix K
% sigma_0: a prior unit weight standard deviation              
Q_mat = (1.0/(sigma_0^2))*K_mat;  % cofactor matrix
P_mat_0 = inv(Q_mat);             % initial weight matrix

disp('Assign the initial weight matrix done');

x_temp=x_0; % assign initial value for unknown vector
P_mat_temp=P_mat_0;
ex_iter_count=0;

outlier_mask=zeros(op_count * scan_count, 1);

% External loop for detecting and removing outliers 
% Introduce the structure temp_adjustment_data
% documentation (comments here)
temp_adjustment_data = struct('x',x_temp, 'y',y, 'P',P_mat_temp, 'sigma_0',...
   sigma_0, 'dt', delta_t, 'op', ops, 'scans', {scans_in_sphe}, 'ap_count',...
   ap_count, 'max_iter_count', max_in_iter, 'incre_ratio_thre', incre_ratio_thre,...
   'outlier_mask', outlier_mask);
 
% solved issue: to put the cell array into a struct, you need to use { }
% outside the cell array to make it a cell
     
% iterate until reaching the termination criteria
% Adjustment calculation
[x_p, Q_xx_mat, res_vec]= RunGMMAdjust(temp_adjustment_data);
 


