function x_0_vec = generate_init(OPs,Scans,s)
% generate_init Calculate initial values for the project
% Unknown parameters x_vec [a0, b1, b2, c0, omega_1, phi_1, kappa_1, Xs_1, Ys_1, Zs_1, ..., omega_s, phi_s, kappa_s, Xs_s, Ys_s, Zs_s]
% Output: initial values for all the unknowns,  'x_0_vec':  4+(s*6) * 1 column vector

%% Initial values for additional calibration parameters (APs)
Ap_0 = [0;0;0;0];

%% Initial values for scanner position and orientations
    
end