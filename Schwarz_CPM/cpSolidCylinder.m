function [cpx,cpy,cpz, dist, bdy] = cpSolidCylinder(x,y,z, zlim, R, cen)
%cpSolidCylinder  Closest point function for a solid cylinder
%   A cylinder rising in the z-direction.
%   [cpx,cpy,cpz, dist, bdy] = cpCylinder(x,y,z, zlim, R, cen)
%     'zlim' defaults to [-1 1]
%     radius 'R' defaults to 1.
%     'cen', location in xy-plane of cylinder, default: [0,0].
%
%   Code is vectorized: any size/shape for x should work.


  % default radius
  if (nargin < 5),   R = 1;   end
  % default bottom/top
  if (nargin < 4),   zlim = [-1  1];   end
  % default center (in x,y) is the origin
  if (nargin < 6),   cen = [0,0];   end

  % shift to the origin
  x = x - cen(1);
  y = y - cen(2);

  zlo = zlim(1);
  zhi = zlim(2);

  bdy1 = (z < zlo);
  bdy2 = (z > zhi);

  [th, r, zp] = cart2pol(x, y, z);
  cpth = th;
  cpr = (r < R).*r + (r >= R).*(R*ones(size(r)));
  cpzp = zp;
  cpzp(bdy1) = zlo;
  cpzp(bdy2) = zhi;

  [cpx, cpy, cpz] = pol2cart(cpth, cpr, cpzp);

  % TODO: color the different bits of the exterior
  bdy = bdy1 | bdy2 | (r >= R);

  dist = sqrt( (x-cpx).^2 + (y-cpy).^2 + (z-cpz).^2 );
  dist(~bdy) = 0;

  % shift back to center
  cpx = cpx + cen(1);
  cpy = cpy + cen(2);
  %cpz = cpz;
