% rbf = tpsK(x,y)
% Computes matrix for thin plate spline kernel K with
% linear polynomials cardinal on (0,0), (1,0), (0,1)
% Calls on: tps
  function rbf = tpsK(x,y)
  % Define points for cardinal polynomials
  ppoints = [0 0; 1 0; 0 1];
  px = [p1(x) p2(x) p3(x)];
  py = [p1(y) p2(y) p3(y)];
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
  for k=1:3
     rbf = rbf + px(:,k)*py(:,k)';
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
