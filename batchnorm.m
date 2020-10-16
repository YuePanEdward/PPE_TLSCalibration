function norm_vec= batchnorm(Mat)
% batchnorm: calculate the norm of a matrix in batch
%   Input:  a matrix Mat (n * m)
%   Output: a vector norm_vec (n * 1), where each element is the norm of a
%   (m * 1) vector

norm_vec= zeros(size(Mat,1),1);

for i=1:size(Mat,1)
    norm_vec(i,1)=norm(Mat(i,:),2); 
end

end