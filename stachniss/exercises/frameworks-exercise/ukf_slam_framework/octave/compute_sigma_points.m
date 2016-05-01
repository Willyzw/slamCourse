function [sigma_points] = compute_sigma_points(mu, sigma)
% Computes the 2n+1 sigma points according to the unscented transform,
% where n is the dimensionality of the mean vector mu.
% The sigma points should form the columns of sigma_points,
% i.e. sigma_points is an nx2n+1 matrix.

global scale;

% Compute lambda
n = length(mu);
num_sig = 2*n+1;
lambda = scale - n;

% TODO: Compute sigma points

sigma_points= ones(n,2*n+1);
sigma_points(:,1)= mu;
sigmasqr= sqrtm(sigma);

for i= 2:n+1
	sigma_points(:,i)= mu+ (sqrt(n+lambda)*sigmasqr(:,i-1));
	sigma_points(:,i+n)= mu- (sqrt(n+lambda)*sigmasqr(:,i-1));
endfor


end
