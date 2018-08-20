function [cpx,cpy,cpz, dist, bdy] = cpBall_right(x,y,z, R, cen,x_left_end)
% cpBall_right Closest point function for a solid ball with left sector cut from
% x_left_end onwards (x_left_end <= cen(1))


%   [cpx,cpy,cpz, dist, bdy] = cpBall_right(x,y,z, R, cen, x_left_end) returns the
%   closest point and distance of (x,y,z) to the solid ball centred at cen
%   with radius R and the left sector being cut from x=x_left_end
%   onwards


  assert(x_left_end <= cen(1));
  
  % defaults
  if (nargin < 4)
    R = 1;
  end
  if (nargin < 5)
    cen = [0, 0, 0];
  end
  
  if (nargin<6)
      x_left_end= -1;
  end
      
%% Calculate the cp points w.r.t the ball with radius R first  

  % Shift to the origin
  x = x - cen(1);
  y = y - cen(2);
  z = z - cen(3);

  [th, phi, r] = cart2sph(x,y,z);
  bdy1 = r >= R; % Temporary bdy
  r(bdy1) = R;
  [cpx,cpy,cpz] = sph2cart(th, phi, r);

  % Shift back
  cpx = cpx + cen(1);
  cpy = cpy + cen(2);
  cpz = cpz + cen(3);
  
 %% Deal with the cut region
 
 % Note that x^2+y^2+z^2=R.^2
 % sectional area (y^2+z^2) from cutting point x=x_left_end 
 
  sec_area=R^2-x_left_end^2; 
  
  I1= cpx <= x_left_end;
  I2= y.^2+z.^2 <=sec_area;  
  
  cpx(I1)=x_left_end;
  cpy(I1 & I2 )= y(I1 & I2 );
  cpz(I1 & I2 )= z(I1 & I2 );
  
  [new_th,new_r]=cart2pol(y(I1 & ~I2),z(I1 & ~I2));
  [cpy(I1 & ~I2), cpz(I1 & ~I2)]=pol2cart(new_th,sqrt(sec_area)); 
  
  %Update the bdy
  bdy=(bdy1 | I1);
  
  %Compute the dist
  dist=sqrt((x-cpx).^2+(y-cpy).^2+(z-cpz).^2);
   
  


  
  