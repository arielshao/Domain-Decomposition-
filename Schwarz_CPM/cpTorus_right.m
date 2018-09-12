function [cpx,cpy,cpz, dist, bdy] = cpTorus_right(x,y,z, R,r,cen,x_left_end)
% cpTorus_right Closest point function for a torus with right sector cut from
% x_left_end onwards (x_left_end <= cen(1))


%   [cpx,cpy,cpz, dist, bdy] = cpTorus_right(x,y,z, R,r, cen, x_left_end) returns the
%   closest point and distance of (x,y,z) to the Torus centred at cen
%   with major radius R and minor radius r. The right sector is cut from x=x_left_end
%   onwards


  assert(x_left_end<=cen(1));

  % defaults
  if (nargin < 4)
    R = 1;
  end
  if (nargin < 5)
    r = 0.4;
  end
  if (nargin < 6)
    cen = [0 0 0];
  end
  


  
  %% Calculate the cp points w.r.t the torus with radius R, r first  
  % shift to the origin
  x = x - cen(1);
  y = y - cen(2);
  z = z - cen(3);

  [cpx,cpy,cpz] =cpTorus(x, y, z,R,r);

  % shift back
  cpx = cpx + cen(1);
  cpy = cpy + cen(2);
  cpz = cpz + cen(3);
  
  
  %% Deal with the cut region

% Note that (R-(x^2+y^2)^1/2)^2+z^2=r^2;
% from cutting point x=x_right_end 
  
     I= cpx <=x_left_end; 
     cpx(I)=x_left_end;
     
   if cpy==0;
      cpy=cpy+1e-4; %deal with the middle point
   else
      cpy=cpy;
   end
     
    [th,phi]=cart2paramTorus(x(I),y(I),z(I),R);
    
    cpy(I)=sign(cpy(I)+eps).*sqrt((R+r*cos(phi)).^2-x_left_end^2);
 
  % Update the bdy
   bdy=I ; 
  
  % Compute the dist
  dist=sqrt((x-cpx).^2+(y-cpy).^2+(z-cpz).^2);