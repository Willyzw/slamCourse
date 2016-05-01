% This is the main unscented Kalman filter SLAM loop. This script calls all the required
% functions in the correct order.
%
% You can disable the plotting or change the number of steps the filter
% runs for to ease the debugging. You should however not change the order
% or calls of any of the other lines, as it might break the framework.
%
% If you are unsure about the input and return values of functions you
% should read their documentation which tells you the expected dimensions.

% Turn off pagination:
close all
clear all
more off;

format long

% Make tools available
addpath('tools');

% Read world data, i.e. landmarks. The true landmark positions are not given to the robot
%landmarks = read_world('../data/world.dat');
load landmarks
% Read sensor readings, i.e. odometry and range-bearing sensor
%data = read_data('../data/sensor_data.dat');
load data
% Initialize belief
mu = zeros(3,1);
sigma = 0.001*eye(3);
map = [];

% For computing lambda
% scale = lambda + dimensionality
global scale;
scale = 3.0;

thetas_pred = [];
thetas_update = [];
x_pred = [];
x_update = [];
y_pred = [];
y_update = [];
eigens_pred = [];
eigens_update = [];

% Perform filter update for each odometry-observation pair read from the
% data file.
for t = 1:size(data.timestep, 2)
    disp('Time step t ='), disp(t)

    % Perform the prediction step of the UKF
%    [mu, sigma, sig_pnts] = prediction_step(mu, sigma, data.timestep(t).odometry);
    [mu, sigma] = prediction_step(mu, sigma, data.timestep(t).odometry);
    sig_pnts = compute_sigma_points(mu, sigma);
    x_pred = [x_pred, mu(1)];
    y_pred = [y_pred, mu(2)];
    thetas_pred = [thetas_pred, mu(3)];
%   eigens_pred = [eigens_pred, min(eig(sigma))];
    eigens_pred = [eigens_pred, (det(sigma))];
    % Perform the correction step of the UKF
    [mu, sigma, map] = raf_correction_step(mu, sigma, data.timestep(t).sensor, map);
    thetas_update = [thetas_update, mu(3)];
    x_update = [x_update, mu(1)];
    y_update = [y_update, mu(2)];
    eigens_update = [eigens_update, (det(sigma))];

    %Generate visualization plots of the current state of the filter

    plot_state(mu, sigma, landmarks, t, map, data.timestep(t).sensor, sig_pnts);
    disp("Current state vector mu ="), disp(mu)
    disp("Map contains the following landmarks:"), disp(map)

sigma(1:3,1:3)

endfor
close all
plot([1:length(thetas_pred)], thetas_pred, 'gx', 'linewidth', 3)
hold on
plot([1:length(thetas_update)], thetas_update, 'ro', 'linewidth', 3)
line([1 330], [-pi -pi])
line([1 330], [-pi/2 -pi/2])
line([1 330], [pi/2 pi/2])
line([1 330], [pi pi])
grid
figure
plot([1:length(x_pred)], x_pred, 'gx', 'linewidth', 3)
hold on
plot([1:length(x_update)], x_update, 'go', 'linewidth', 3)
plot([1:length(y_pred)], y_pred, 'rx', 'linewidth', 3)
plot([1:length(y_update)], y_update, 'ro', 'linewidth', 3)
grid
%{
figure
hold on
plot([1:length(eigens_pred)], eigens_pred, 'gx', 'linewidth', 3)
plot([1:length(eigens_update)], eigens_update, 'ro', 'linewidth', 3)
%}

disp("Final system covariance matrix:"), disp(sigma)
% Display the final state estimate
disp("Final robot pose:")
disp("mu_robot = "), disp(mu(1:3)), disp("sigma_robot = "), disp(sigma(1:3,1:3))
