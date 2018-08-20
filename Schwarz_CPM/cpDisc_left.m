function [cpx, cpy, dist,bdy] = cpDisc_left(x, y, R, cen, x_right_end)
%cpDisc_left Closest point function for a disc with right sector cut from
%x=x_right_end onwards (x_right_end >= cen(1))


%   [cpx,cpy, dist, bdy] = cpDisc_left(x,y, R, cen, x_right_end) returns the
%   closest point and distance of (x,y) to the disc centred at cen
%   with radius R and the right sector being cut from x=x_right_end
%   onwards

 assert(x_right_end >= cen(1));
 
 % Defaults
  if (nargin < 3)
    R = 1;
  end
  if (nargin < 4)
    cen = [0 0];
  end
  if (nargin<5)
      x_right_end=1;
  end

 %% Calculate the cp points w.r.t the disc with radius R first
  
  % Shift to the origin
  x=x-cen(1);
  y=y-cen(2);
  
  [th,r]=cart2pol(x,y);
  [cpx,cpy]=pol2cart(th,R);
  dist =r-R;
  I= dist <0;
  
  cpx(I)=x(I);
  cpy(I)=y(I);
  dist(I)=0;
  bdy1= ~I; 
  
  % Shift back 
  cpx=cpx+cen(1);
  cpy=cpy+cen(2);
  
  %% Deal with the cut region
  
 % Compute the upper bound for y on the cutting line x=x_right_end
  y_right_end=sqrt(R^2-x_right_end^2);
  
  I1= cpx >x_right_end;
  I2= y> -y_right_end;
  I3= y< y_right_end;
  I4= y<-y_right_end;
  I5= y> y_right_end;
  
  cpx(I1)=x_right_end;
  cpy(I1 & I2 & I3)=y(I1 & I2 & I3);
  cpy(I1 & I4)=-y_right_end;
  cpy(I1 & I5)=y_right_end;
  
 % Update the bdy
   bdy=(I1 | bdy1);
   
 % Compute the dist
   dist=sqrt((x-cpx).^2+(y-cpy).^2);

