% Computes the total error of the graph
function Fx = compute_global_error(g)

Fx = 0;

% Loop over all edges
for eid = 1:length(g.edges)
  edge = g.edges(eid);

  % pose-pose constraint
  if (strcmp(edge.type, 'P') ~= 0)

    x1 = g.x(edge.fromIdx:edge.fromIdx+2);  % the first robot pose
    x2 = g.x(edge.toIdx:edge.toIdx+2);      % the second robot pose
    X1 = v2t(x1);
    X2 = v2t(x2);
    
    %TODO compute the error of the constraint and add it to Fx.
    % Use edge.measurement and edge.information to access the
    % measurement and the information matrix respectively.
    z12 = edge.measurement;
    Z12 = v2t(z12);
%     error = Z12(1:2,1:2)' * (X1(1:2,1:2)' * (x2(1:2)-x1(1:2)) - z12(1:2));
%     error(3) = x2(3) - x1(3) - z12(3)
    error = t2v(inv(Z12)*(inv(X1)*X2));
    Fx = Fx + error'*edge.information*error;

  % pose-landmark constraint
  elseif (strcmp(edge.type, 'L') ~= 0)
    x = g.x(edge.fromIdx:edge.fromIdx+2);  % the robot pose
    l = g.x(edge.toIdx:edge.toIdx+1);      % the landmark

    %TODO compute the error of the constraint and add it to Fx.
    % Use edge.measurement and edge.information to access the
    % measurement and the information matrix respectively.
    X = v2t(x);
    z = edge.measurement;
    error = X(1:2,1:2)'*(l-x(1:2)) - z;
    Fx = Fx + error'*edge.information*error;
  end

end
