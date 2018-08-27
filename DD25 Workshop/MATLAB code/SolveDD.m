function [u,err]=SolveDD(f,eta,a,b,gg,gd,Ii,d,u0,CC,N)
% SOLVEDD solves 1d model problem with domain decomposition
%   u=SolveDD(f,eta,a,b,gg,gd,Ii,d,u0,CC,N); solves (eta-dxx)u=f on
%   the 1d domain (a,b) with Dirichlet conditions u=gg and u=gd on the
%   equidistant grid defined by the length of the rhs f using a
%   Schwarz algorithm. Subdomains are defined by the indices Ii(j) of
%   a non-overlapping decomposition, enlarged by d mesh sizes in both
%   directions to obtain an overlapping one, using the initial guess
%   u0, doing N iterations. CC==0 uses no coarse grid, CC==1 a coarse
%   grid with nodes centered in the subdomains, and CC==2 an optimized
%   coarse grid for parallel Schwarz or RAS

ue=Solve1d(f,eta,a,b,gg,gd);                % reference solution
u=[gg;u0;gd];                               % initial guess
J=length(f); h=(b-a)/(J+1); x=a:h:b;
ai=[1,Ii(2:end-1)-d]; bi=[Ii(2:end-1)+d,J]; % construct overlapping dec
A=A1d(eta,a,b,J);
err(1)=norm(u-ue,inf);
ft=f; ft(1)=ft(1)+1/h^2*gg; ft(end)=ft(end)+1/h^2*gd;
if CC==1                                    % compute coarse components
  Im=Coarse(Ii); [E,R,Ac]=CoarseOperators(Im,A); 
elseif CC==2
  Im=CoarseOpt(Ii); [E,R,Ac]=CoarseOperators(Im,A);  
end;

for n=1:N                                   % Schwarz iterations
  uo=u;
  for j=1:length(ai)                        % subdomain solves
    tmp=Solve1d(f(ai(j):bi(j)),eta,(ai(j)-1)*h,(bi(j)+1)*h,uo(ai(j)),uo(bi(j)+2)); 
    if j==1                                 % compose a global solution
      u(ai(j):bi(j)+2-d+1)=tmp(1:end-d+1);
    elseif j==length(ai)
      u(ai(j)+d+1:bi(j)+2)=tmp(1+d+1:end);
    else  
      u(ai(j)+d+1:bi(j)+2-d+1)=tmp(1+d+1:end-d+1);
    end
  end;
  if CC                                     % coarse grid correction
    r=ft-A*u(2:end-1);                      % compute residual
    uc=Ac\(R*r);                            % compute coarse correction
    u(2:end-1)=u(2:end-1)+E*uc;             % add coarse correction 
  end; 
  err(n+1)=norm(u-ue,inf);      
  plot(x,u,'o',x,ue,'-');                   % plot solution and approximation
  title(['iteration number ' num2str(n)])
  xlabel('x')
  pause
end;
