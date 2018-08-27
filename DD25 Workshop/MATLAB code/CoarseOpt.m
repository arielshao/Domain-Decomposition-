function Im=CoarseOpt(Ii,d)
% COARSEOPT compute optimized coarse grid location in one dimension
%   Im=Coarse(Ii,d) computes for a non overlapping domain decomposition
%   given by the vector of interface indices Ii coarse grid nodes
%   located around the discontinuities of a parallel Schwarz method
%   if d is given, the coarse nodes are placed close to the
%   interfaces, so one can see the correction in the overlap is not
%   optimal. Otherwise the coarse correction is perfectly placed
%   next to the non-overlapping interfaces and gives the solution
%   directly. 

if nargin<2
  for i=1:length(Ii)-2;                   % coarse grid points in
    Im(2*i-1)=Ii(i+1)-1;                  % the center of overlaps
    Im(2*i)=Ii(i+1);                      % permitting discontinuities 
  end;
else 
  for i=1:length(Ii)-2;                   % coarse grid points next
    Im(2*i-1)=Ii(i+1)-d;                  % to the interfaces, still
    Im(2*i)=Ii(i+1)+d-1;                  % optimal but not in the overlap 
  end;
end
