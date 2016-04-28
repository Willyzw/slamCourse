more off;
clc;
clear all;
close all;
format long g;

addpath('tools');

% load the graph into the variable g
% only leave one line uncommented

% simulation datasets
load ../data/simulation-pose-pose.mat
% load ../data/simulation-pose-landmark.mat

% real-world datasets
% load ../data/intel.mat
% load ../data/dlr.mat
% g.x = g.x(1:12);
% g.edges = g.edges(1:3);
% plot the initial state of the graph
% plot_graph(g, 0);

fprintf('Initial error %f\n', compute_global_error(g));
% compute_global_error(g)

% the number of iterations
numIterations = 100;

% maximum allowed dx
EPSILON = 10^-4;

% Error
err = 0;

% carry out the iterations
for i = 1:numIterations
  fprintf('Performing iteration %d\n', i);

  dx = linearize_and_solve(g);

  % TODO: apply the solution to the state vector g.x
  g.x(4:end) = g.x(4:end) + dx(4:end);
  
  % plot the current state of the graph
  plot_graph(g, i);

  err = compute_global_error(g);

  % Print current error
  fprintf('Current error %f\n', err);

  % TODO: implement termination criterion as suggested on the sheet
  if (max(abs(dx)) < EPSILON) 
      break;
  end
end

fprintf('Final error %f\n', err);
