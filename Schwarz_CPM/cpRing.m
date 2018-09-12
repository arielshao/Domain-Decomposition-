function [cpx, cpy, dist, bdy ] = cpRing(x,y,r,R,cen);

% cpRing is the closest point function a ring with inner radius r, outer
% radies R;
if (nargin < 3)
   r=0.5;
end
if (nargin < 4)
   R=1;
end
if (nargin < 5)
    cen = [0 0];
end

  % shift to the origin
  x = x - cen(1);
  y = y - cen(2);
 [th, rr] = cart2pol(x, y);
 [cpx1, cpy1] = pol2cart(th, R);
 [cpx2, cpy2] = pol2cart(th, r);
 dist_R = rr - R;  %negative inside disc with big radius R
 dist_r = rr -r;   %negative inside disc with small radius r
 I1=dist_R>0;
 I2=dist_r <0;
 I = (dist_R<0) & (dist_r>0);
 cpx=cpx1;
 cpy=cpy1;
 cpx(I2)=cpx2(I2);
 cpy(I2)=cpy2(I2);
 cpx(I) = x(I);
 cpy(I) = y(I);
 bdy=~I;
 dist(I)=0;
 dist(I1)=dist_R(I1);
 dist(I2)=-dist_r(I2);
 
  % shift back
  cpx = cpx + cen(1);
  cpy = cpy + cen(2);







