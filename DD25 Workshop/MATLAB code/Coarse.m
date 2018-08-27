function Im=Coarse(Ii)
% COARSE compute coarse grid location in one dimension
%   Im=Coarse(Ii) computes for a non overlapping domain decomposition
%   given by the vector of interface indices Ii coarse grid nodes
%   located in the center of subdomains

I=length(Ii)-1;                         
for i=1:I
  Im(i)=round((Ii(i+1)+Ii(i))/2);
end;
