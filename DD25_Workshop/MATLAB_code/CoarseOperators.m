function [E,R,Ac]=CoarseOperators(Im,A)
% COARSE compute coarse grid components in one dimension
%   [E,R,Ac]=CoarseOperators(Im,A), computes for a given coarse grid
%   Im and the discretization matrix A a linear extension operator E,
%   the corresponding restriction operator R, and the coarse system
%   matrix Ac.

J=size(A,1); I=length(Im);
R=sparse(I-1,J); E=sparse(J,I-1);
for i=1:I
  if i==1
    dx=1/Im(i); E(1:Im(i)-1,i)=dx:dx:1-dx; % linear extension
  else
    dx=1/(Im(i)-Im(i-1)); E(Im(i-1):Im(i)-1,i)=0:dx:1-dx;  
  end
  if i==I
    dx=1/(J-Im(i)+1); E(Im(i):J,i)=1:-dx:dx;         
  else
    dx=1/(Im(i+1)-Im(i)); E(Im(i):Im(i+1)-1,i)=1:-dx:dx;         
  end
end;
R=E'; R=spdiags(1./sum(R')',0,I,I)*R;      % restriction
Ac=R*A*E;                                  % Galerkin coarse matrix
