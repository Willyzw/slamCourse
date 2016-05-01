function particles = prediction_step(particles, u, noise)
% Updates the particles by drawing from the motion model
% Use u.r1, u.t, and u.r2 to access the rotation and translation values
% which have to be pertubated with Gaussian noise.
% The position of the i-th particle is given by the 3D vector
% particles(i).pose which represents (x, y, theta).

% noise parameters
% Assume Gaussian noise in each of the three parameters of the motion model.
% These three parameters may be used as standard deviations for sampling.
r1Noise = noise(1);
transNoise = noise(2);
r2Noise = noise(3);

numParticles = length(particles);

for i = 1:numParticles

  % append the old position to the history of the particle
  particles(i).history{end+1} = particles(i).pose;

  % TODO: sample a new pose for the particle
  r1 = normrnd(u.r1, r1Noise);
  r2 = normrnd(u.r2, r2Noise);
  trans = normrnd(u.t, transNoise);
  oldPose = particles(i).pose;
  particles(i).pose(1) = oldPose(1) + trans*cos(oldPose(3) + r1);
  particles(i).pose(2) = oldPose(2) + trans*sin(oldPose(3) + r1);
  particles(i).pose(3) = oldPose(3) + r1 + r2;

end

end
