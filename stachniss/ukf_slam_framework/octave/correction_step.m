%Nichola
function [mu, sigma, map] = correction_step(mu, sigma, z, map);
% Updates the belief, i.e., mu and sigma after observing landmarks,
% and augments the map with newly observed landmarks.
% The employed sensor model measures the range and bearing of a landmark
% mu: state vector containing robot pose and poses of landmarks obeserved so far.
% Current robot pose = mu(1:3)
% Note that the landmark poses in mu are stacked in the order by which they were observed
% sigma: the covariance matrix of the system.
% z: struct array containing the landmark observations.
% Each observation z(i) has an id z(i).id, a range z(i).range, and a bearing z(i).bearing
% The vector 'map' contains the ids of all landmarks observed so far by the robot in the order
% by which they were observed, NOT in ascending id order.

% For computing sigma
global scale;

% Number of measurements in this time step
m = size(z, 2);

% Measurement noise
Q = 0.01*eye(2);

for i = 1:m
	% Get the id of the landmark corresponding to the i-th observation
	landmarkId = z(i).id;

	% Get the actual measurement made
	z_actual = [z(i).range; z(i).bearing];

	% If the landmark is observed for the first time:
	if(isempty(find(map==landmarkId)))
	    % Add it to the map
	    map = [map; landmarkId];
	    % TODO: Initialize its pose according to the measurement and add it to mu
%	    newX = mu(1) + z(i).range*cos(mu(3,:) + z(i).bearing);
%	    newY = mu(2) + z(i).range*sin(mu(3,:) + z(i).bearing);
%	    mu = [mu; newX; newY];
	    mu = [mu; z(i).range(); z(i).bearing()];
	    % TODO: Initialize its uncertainty and add it to sigma
%	    sigma = blkdiag(sigma, 1*eye(2));
	    sigma = blkdiag(sigma, Q);

	    sig_pnts_new = compute_sigma_points(mu, sigma);
	    sig_pnts_new(3,:) = normalize_angle(sig_pnts_new(3,:));
	    newX = sig_pnts_new(1,:) + sig_pnts_new(end-1,:).*cos(sig_pnts_new(3,:) + sig_pnts_new(end,:));
	    newY = sig_pnts_new(2,:) + sig_pnts_new(end-1,:).*sin(sig_pnts_new(3,:) + sig_pnts_new(end,:));
	    sig_pnts_new(end-1,:) = newX;
	    sig_pnts_new(end,:) = newY;

	    n = length(mu);
	    lambda = scale - n;
	    w0 = lambda/scale;
	    wm = [w0, repmat(1/(2*scale),1,2*n)];
	    cosines = sum(cos(sig_pnts_new(3,:)).*wm);
	    sines = sum(sin(sig_pnts_new(3,:)).*wm);
	    mu_theta = normalize_angle(atan2(sines,cosines));
	    idx = 2:size(sig_pnts_new,2);
	    mu = (2*lambda*sig_pnts_new(:,1) + sum(sig_pnts_new(:,idx),2)) / (2*scale); 
	    mu(3) = mu_theta;
	
	    diff = sig_pnts_new - repmat(mu,1,size(sig_pnts_new,2));
	    diff(3,:) = normalize_angle(diff(3,:));
	    sigma = (2*lambda*diff(:,1)*diff(:,1)' + diff(:,idx)*diff(:,idx)') / (2*scale);
	    continue;
	endif

	% TODO: Compute sigma points from the predicted mean and covariance
	sigma_points = compute_sigma_points(mu, sigma);
	sigma_points(3,:) = normalize_angle(sigma_points(3,:));

	% Compute lambda
	n = length(mu);
	num_sig = size(sigma_points,2);
	lambda = scale - n;

	% TODO: Compute z_points (2x2n+1), which consists of predicted measurements from all sigma points
	loc = find(map==(z(i).id));
	landmarkXs = sigma_points(2*loc + 2,:);
	landmarkYs = sigma_points(2*loc + 3,:);
	z_points_range = sqrt((landmarkXs - sigma_points(1,:)).^2 + (landmarkYs - sigma_points(2,:)).^2);
	z_points_bearing = atan2(landmarkYs-sigma_points(2,:), landmarkXs-sigma_points(1,:));
	z_points_bearing = z_points_bearing - sigma_points(3,:);
	z_points_bearing = normalize_angle(z_points_bearing);
	z_points = [z_points_range;z_points_bearing];

	% TODO: Compute the weights of the sigma points and compute zm.
	% zm is the recovered expected measurement mean from z_points.
	% It will be a 2x1 vector [expected_range; expected_bearing].
	zm_range = (lambda*z_points_range(1) + 0.5*sum(z_points_range(2:end), 2)) / scale;
	w0 = lambda/scale;
	wm = [w0, repmat(1/(2*scale),1,2*n)];

	cosines = sum(cos(z_points_bearing).*wm);
	sines = sum(sin(z_points_bearing).*wm);
	zm_bearing = normalize_angle(atan2(sines,cosines));
	zm = [zm_range; zm_bearing];

	% TODO: Compute the innovation covariance matrix S (2x2).
	dz = z_points - repmat(zm,1,num_sig);
	dz(2,:) = normalize_angle(dz(2,:));

	Pzz = (2*lambda*dz(:,1)*dz(:,1)' + dz(:,2:end)*dz(:,2:end)') / (2*scale);
	S = Pzz + Q;
%	Sc = chol(S);  % note: S = Sc'*Sc
%	Sci = inv(Sc);  % note: inv(S) = Sci*Sci'

	% TODO: Compute sigma_x_z (which is equivalent to sigma times the Jacobian H transposed in EKF).
	% sigma_x_z is an nx2 matrix, where n is the current dimensionality of mu
	dx = sigma_points - repmat(mu,1,num_sig);
	dx(3,:) = normalize_angle(dx(3,:));
	sigma_x_z = (2*lambda*dx(:,1)*dz(:,1)' + dx(:,2:end)*dz(:,2:end)') / (2*scale);

	% TODO: Compute the Kalman gain
%	Wc = sigma_x_z * Sci;
%	K = Wc * Sci';
	K = sigma_x_z*inv(S);

	% TODO: Update mu and sigma
	diff_z = (z_actual - zm);
	diff_z(2,:) = normalize_angle(diff_z(2,:));
	
	mu = mu + K*diff_z;
%	sigma = sigma - Wc*Wc';
	sigma = sigma - K*S*K';
	% TODO: Normalize the robot heading mu(3)
	mu(3) = normalize_angle(mu(3));

%	disp('sigma new after'), disp(sigma(end-1:end,end-1:end));
	%PSD_check = chol(sigma);
endfor
end
