function [cpx, cpy, dist,bdy] = cpRectangleDomain(x, y, a, b, cen)
%cpRectangleDomain Closest point function for a rectangle domain
%   [cpx, cpy, dist] = cpRectangleDomain(x, y)
%      A rectangle domain with horizontal length a, and vertical length b centered at
%      the origin.
%   [cpx, cpy, sdist] = cpRectangleDomain(x, y, xc, yc)
%       A rectangle with horizontal side length a, and vertical length b centered at
%      [xc, yc]
%   Not vectorized, uses a loop internally.

  % defaults
  if (nargin < 5)
    cen = [0, 0];
  end
  
  if (nargin<3)
     a=2;b=1;
  end

  % shift to the origin
  x=x - cen(1);
  y=y- cen(2);

  cpx = zeros(size(x));
  cpy = zeros(size(y));

  for i=1:length(x(:))
    xx = x(i);
    yy= y(i);

    if (yy >= b/2)
      cpyy = b./2;
    elseif (yy <= -b/2)
      cpyy = -b/2;
    else
      cpyy = yy;
    end

    if (xx >= a/2)
      cpxx = a/2;
    elseif (xx <= -a/2)
      cpxx = -a/2;
    else
      cpxx = xx;
    end
    cpx(i) = cpxx;
    cpy(i) = cpyy;
  end
  
   % Compute the dist
  dist=sqrt((x-cpx).^2+(y-cpy).^2);
  
 % points outside the rectangle domain
  I1=(dist>0); 
 % points on the rectangle boundary
  I2=(cpx==a/2 | cpx== -a/2);
  I3=(cpy==b/2 | cpy== -b/2);
  
 %Update the bdy
  bdy=(I1 | I2 | I3);
  
 % shift back
  cpx = cpx + cen(1);
  cpy = cpy + cen(2);
 
  
  
  

