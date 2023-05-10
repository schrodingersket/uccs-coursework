% tf = testfunction(s,points)
% Evaluates testfunction
% prod_{d=1}^s x_d*(1-x_d)  (normalized so that its max is 1)
% at s-dimensional points
function tf = testfunction(s,points)
tf = 4^s*prod(points.*(1-points),2);
