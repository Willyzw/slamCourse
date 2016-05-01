
mu = [2;1]
sigma = [0.1 0.0; 0.0 0.1]


function [res] = g(x)

%res = x + 1;
res(1) = 1 + x(1) + sin(2*x(1)) * cos(x(2));
res(2) = 0.2*x(2) + 2.0; 
end


alpha = 0.25;
kappa = 10;
beta = 2;

axis([-3 3 -3 3],"square");
hold off

n = length(mu);


lambda = alpha*alpha*(n+kappa)-n;

[sig_points, w_m, w_s] = compute_sigma_points(mu, sigma, lambda, alpha, beta);


plot(mu(1),mu(2),'ro','markersize',8)
hold on

drawprobellipse(mu, sigma, 0.9, 'r');
plot(sig_points(1,:),sig_points(2,:),'rx','markersize',8)


for i=1:length(sig_points)
    sig_points(:,i) = g(sig_points(:,i));
end


plot(sig_points(1,:),sig_points(2,:),'kx','markersize',8)


w = zeros(size(sig_points));
w(:,1) = w_m(1);
w(:,2:end) = w_m(2);
temp = sig_points.*w;
mu = sum(temp,2);
plot(mu(1),mu(2),'ko','markersize',8);

diff = sig_points-repmat(mu,1,2*n+1);

sigma = zeros(n);
for i=1:2*n+1
sigma = sigma + w_s(i)*diff(:,i)*diff(:,i)';
endfor
drawprobellipse(mu, sigma, 0.9, 'k');

mu
sigma