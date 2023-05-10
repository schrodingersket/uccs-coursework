% rbf = tpsH(x,y,a)
% Computes matrix for homogeneous thin plate spline kernel kappa with
% linear polynomials cardinal on (0,0), (a,0), (0,a)
% Inputs x and y are assumed to be scaled to [0,a]^2
% Calls on: tps
  function rbf = tpsH(x,y,a)
  % Define points for cardinal polynomials
  ppoints = a*[0 0; 1 0; 0 1];
  px = [p1(x/a) p2(x/a) p3(x/a)];
  py = [p1(y/a) p2(y/a) p3(y/a)];
  r = DistanceMatrix(x,y); rbf = tps(1,r);
  for k=1:3
     r = DistanceMatrix(ppoints(k,:),y);
     rbf = rbf - px(:,k)*tps(1,r);
     r = DistanceMatrix(x,ppoints(k,:));
     rbf = rbf - tps(1,r)*py(:,k)';
  end
  for j=1:3
     for k=1:3
        r = DistanceMatrix(ppoints(j,:),ppoints(k,:));
        rbf = rbf + px(:,j)*py(:,k)'*tps(1,r);
     end
  end
  return
  % The cardinal polynomials
  function w = p1(z)
  w = 1 - z(:,1) - z(:,2);
  return
  function w = p2(z)
  w = z(:,1);
  return
  function w = p3(z)
  w = z(:,2);
  return
