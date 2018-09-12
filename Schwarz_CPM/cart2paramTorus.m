function [az,elev,r] = cart2paramTorus(x,y,z,R)
%CART2PARAMTORUS Transform Cartesian to the spherical parametrization of a torus
% with major radius R;
% [TH,PHI] = CART2PARAMTORUS(X,Y,Z,R) transforms corresponding elements of
% data stored in Cartesian coordinates X,Y,Z to spherical parametrization
% of a torus (azimuth TH, elevation PHI).  The arrays
% X,Y, and Z must be the same size (or any of them can be scalar).
% TH and PHI are returned in radians.
%
% TH is the counterclockwise angle in the xy plane measured from the
% positive x axis.  PHI is the elevation angle from the xy plane.
hypotxy = hypot(x,y);
r= hypot(z,hypotxy);
elev = atan2(z,hypotxy-R);
az = atan2(y,x);
