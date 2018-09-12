function [cpx, cpy, dist,bdy] = cpEllipticalDisc_right(x, y,a,b, cen, x_left_end)
%cpEllipticalDisc_right Closest point function for a elliptical disc with 
%left sector cut from x=x_left_end onwards (x_left_end >= cen(1))


%   [cpx,cpy, dist, bdy] = cpEllipticalDisc_right(x,y, a,b, cen, x_left_end) returns the
%   closest point and distance of (x,y) to the elliptical-shape disc centred at cen
%   with major axis 'a' and minor axis 'b'.  and the left sector being cut 
%   from x=x_left_end onwards

 assert(x_left_end <= cen(1));
 
% Defaults
  if (nargin < 4)
    a=1;b=1
  end
  if (nargin < 5)
    cen = [0 0];
  end
  if (nargin<6)
      x_left_end=-a;
  end

 %% Find the closest points to the elliptical-shaped disc first
  [cpx,cpy,dist,bdy1]=cpEllipticalDisc(x,y,a,b,cen);

 %% Deal with the cutted region
 
  y_left_end=b*sqrt(1-x_left_end^2/a^2);
  I1= cpx<=x_left_end;
  I2= y>= -y_left_end;
  I3= y<= y_left_end;
  I4= y<-y_left_end;
  I5= y> y_left_end;
  
  cpx(I1)=x_left_end;
  cpy(I1 & I2 & I3)=y(I1 & I2 & I3);
  cpy(I1 & I4)=-y_left_end;
  cpy(I1 & I5)=y_left_end;

% Update the bdy
   bdy=(I1 | bdy1);
 % Compute the dist
dist=sqrt((x-cpx).^2+(y-cpy).^2);
