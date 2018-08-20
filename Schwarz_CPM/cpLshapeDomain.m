function [cpx, cpy, dist,bdy] = cpLshapeDomain(x, y, a, b, cen)
%cpLshapeDomain  Closest point function for a Lshape domain which consists of two
% identical rectangles with longer side of length a, and shorter side of
% length b
%   [cpx, cpy, dist, bdy] = cpLshapeDomain(x, y)
%   A L-shape domain with the lower rectangle of horizontal length 2, 
%   and vertical length 1 centered at the origin.
%   [cpx, cpy, dist,bdy] = cpLshapeDomain(x, y, xc, yc)
%   A L-shape domain with the lower rectangle of horizontallength 2, 
%   and vertical length 1 centered at [xc, yc]
%
%   Note: returns signed distance (with negative inside).
%
%   Not vectorized, uses a loop internally.
  
 assert(a>b);
  % defaults
  if (nargin < 5)
    cen = [0, 0];
  end

  if (nargin<3)
     a=2;b=1;
  end
  
  % shift to the origin
  x = x - cen(1);
  y = y - cen(2);

  cpx = zeros(size(x));
  cpy = zeros(size(y));

  for i=1:length(x(:))
    xx = x(i);
    yy = y(i);

    if (yy >= a-b & xx>=0)
      cpyy = a-b;
    elseif (yy>=a-b & xx<=0 & yy>=-xx)
      cpyy= a-b;
    elseif (yy <= -b)
      cpyy = -b;
    elseif(yy>=0 & yy<-xx)
      cpyy=0;
    else
      cpyy = yy;
    end

    if (xx >= b)
      cpxx = b;
    elseif (xx <= b-a & yy<=0 )
      cpxx = b-a;
    elseif (xx<= b-a & yy>0 & yy<-xx)
        cpxx= b-a;
    elseif (xx<=0 & yy>=0 & yy>=-xx)
        cpxx=0;
     else
      cpxx = xx;
    end
    cpx(i) = cpxx;
    cpy(i) = cpyy;
  end
  
% Euclidean distance between the original points and their closest points
  dist = sqrt( (x-cpx).^2 + (y-cpy).^2 );
 
   % points outside the Lshape domain
  I1=(dist>0); 
 % points on the rectangle boundary
  I2=(cpx==b | cpx== b-a );
  I3=( cpx>b-a & cpx<=0 & cpy==0);
  I4=(cpx==0 &  cpy>=0 & cpy<=a-b);
  I5=(cpy==-b |cpy==a-b);
  
  bdy=(I1 | I2 | I3| I4 |I5);
  % shift back
  cpxx = cpxx + cen(1);
  cpyy = cpyy + cen(2);

