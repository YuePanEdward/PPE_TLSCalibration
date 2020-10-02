function P_mat = stochastic_model(sigma_0,sigma_rho,sigma_theta,sigma_alpha,n)
% stochastic_model£ºDefine the stochastic model and output the weight matrix
% Input: a prior standard deviation: sigma_0
%        measurement standard deviation for range, horizotal angle and vertical angle: sigma_rho,sigma_theta,sigma_alpha
%        total number of measurments: n (scan number times op number)
% Output: P_mat: the weight matrix


D_mat = diag(repmat([sigma_rho,sigma_theta,sigma_alpha].^2,1,n)); % vcv-matrix
Q_mat = D_mat./(sigma_0^2); % cofactor matrix
P_mat = inv(Q_mat);         % weight matrix 

end