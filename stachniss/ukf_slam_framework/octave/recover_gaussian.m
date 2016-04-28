function [mu, sigma] = recover_gaussian(sigma_points, w_m, w_c)

w = zeros(size(sigma_points));
w(:,1) = w_m(1);
w(:,2:end) = w_m(2);
temp = sigma_points.*w;
mu = sum(temp,2);

n = length(mu);

diff = sigma_points-repmat(mu,1,2*n+1);

sigma = zeros(n);
for i=1:2*n+1
sigma = sigma + w_c(i)*diff(:,i)*diff(:,i)';
endfor

