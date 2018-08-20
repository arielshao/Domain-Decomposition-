function [cpx,cpy,cpz, dist, bdy] = cpBall(x,y,z, R, cen)
%cpBall  Closest point function for a solid ball.
%   [cpx,cpy,cpz, dist, bdy] = cpBall(x,y,z, R) returns the
%   closest point and distance to (x,y,z).  If R is omitted it
%   defaults to a unit ball, centered at the origin.
%
%   [cpx,cpy,cpz, distm bdy] = cpBall(x,y,z, R, [xc,yc,zc]) is a ball
%   of radius R, centered at (xc,yc,zc)
%

  % defaults
  if (nargin < 4)
    R = 1;
  end
  if (nargin < 5)
    cen = [0, 0, 0];
  end

  % shift to the origin
  x = x - cen(1);
  y = y - cen(2);
  z = z - cen(3);

  [th, phi, r] = cart2sph(x,y,z);
  bdy = r >= R;
  dist = r - R;
  r(bdy) = R;
  [cpx,cpy,cpz] = sph2cart(th, phi, r);

  %dist = sqrt( (x-cpx).^2 + (y-cpy).^2 + (z-cpz).^2 );
  dist(~bdy) = 0;

  % shift back
  cpx = cpx + cen(1);
  cpy = cpy + cen(2);
  cpz = cpz + cen(3);
