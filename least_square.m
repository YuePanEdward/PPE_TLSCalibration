function unknown_diff  = least_square( Amat, bvec, Pmat )
%Get the change of unknown parameters by least square
%   Detailed explanation goes here

unknown_diff = (Amat' * Pmat * Amat)^(-1) * Amat * Pmat * bvec;

end

