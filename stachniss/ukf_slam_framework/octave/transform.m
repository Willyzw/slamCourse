function [points] = transform(points)


% Linear function
points(1,:) = points(1,:) + 1;
points(2,:) = points(2,:) + 2;


%{
% Nonlinear function 1
x = points(1,:);
y = points(2,:);
r = sqrt(sum([x.*x; y.*y]));
theta = atan2(y,x);
points = [r;theta];
%}


%{
% Nonlinear function 2
points(1,:) = points(1,:).*cos(points(1,:)).*sin(points(1,:));
points(2,:) = points(2,:).*cos(points(2,:)).*sin(points(2,:));
%}
