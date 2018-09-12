function [cpx,cpy,cpz, dist, bdy] = cpSphere_left(x,y,z, R, cen,x_right_end)
% cpSphere_left Closest point function for a sphere with right sector cut from
% x_right_end onwards (x_right_end >= cen(1))


%   [cpx,cpy,cpz, dist, bdy] = cpSphere_left(x,y,z, R, cen, x_right_end) returns the
%   closest point and distance of (x,y,z) to the sphere centred at cen
%   with radius R and the right sector being cut from x=x_right_end
%   onwards


  assert(x_right_end >= cen(1));

  % defaults
  if (nargin < 4)
    R = 1;
  end
  if (nargin < 5)
    cen = [0, 0, 0];
  end

  
  %% Calculate the cp points w.r.t the sphere with radius R first  
  % shift to the origin
  x = x - cen(1);
  y = y - cen(2);
  z = z - cen(3);

  [th, phi, r] = cart2sph(x,y,z);
  [cpx,cpy,cpz] = sph2cart(th, phi, R);


  % shift back
  cpx = cpx + cen(1);
  cpy = cpy + cen(2);
  cpz = cpz + cen(3);
  
  
  %% Deal with the cut region

% Note that x^2+y^2+z^2=R.^2
% sectional area (y^2+z^2) from cutting point x=x_right_end 

  sec_area=R^2-x_right_end^2; 
  
  I= cpx >=x_right_end; 
  
  cpx(I)=x_right_end;
  
  [new_th,new_r]=cart2pol(y(I),z(I));
  [cpy(I), cpz(I)]=pol2cart(new_th,sqrt(sec_area));
  
  % Update the bdy
  bdy=I ; 
  
  % Compute the dist
  dist=sqrt((x-cpx).^2+(y-cpy).^2+(z-cpz).^2);