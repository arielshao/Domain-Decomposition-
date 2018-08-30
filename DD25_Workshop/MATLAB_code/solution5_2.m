RoomData                                    % include problem parameters
a=9; %d=11;                                  % a and a+d are the interfaces
d=12;
maxiter=20; 
f1=f(:,1:a+d-1); f2=f(:,a+1:end);           % subdomain source terms
z=zeros(J,1);
G1=zeros(a+d+1,1);G1(end-d)=1;              % trace operators
G2=zeros(J-a+2,1);G2(d+1)=1;                
b=[Solve2d(f1,eta,0,a+d,gg,z)*G1;           % construct interface system    
   Solve2d(f2,eta,a,J+1,z,gd)*G2];          % T*g=b, g unknowns at interfaces
T=@(g) [g(1:J)-Solve2d(zeros(size(f1,1),size(f1,2)),eta,0,a+d,z,g(J+1:2*J))*G1;
        g(J+1:2*J)-Solve2d(zeros(size(f2,1),size(f2,2)),eta,a,J+1,g(1:J),z)*G2];
g=zeros(2*J,1);                             % zero initial guess
ri(1)=norm(b-T(g));                         % initial residual
for i=1:maxiter                             
  g=g-T(g)+b;                               % Schwarz iteration
  ri(i+1)=norm(b-T(g));                     % keep residual for plotting
end
[g,fl,r,it,rk]=gmres(T,b);                  % use GMRES Krylov solver
semilogy(0:length(ri)-1,ri,'-o',0:length(rk)-1,rk,'-+');
xlabel('iteration')
ylabel('residual')
legend('Iterative','Krylov')

%% Observations and Analysis 

%       d     Krylov iterations (error)     Schwarz Iterations(error)
%   -----------------------------------------------------------------------------
%       2         10       1.23e-6           20  1.13e-3      
%       3         8        1.05e-6           20  5.85e-5
%       4         8        1.37e-7           20  2.60e-6       
%       5         7        4.27e-7           7   4.50e-3
%       6         7        3.21e-8           7   1.40e-3
%       7         6        3.29e-7           6   1.39e-3
%       8         6        4.66e-8           6   4.28e-4
%       9         5        4.95e-7           5   4.51e-4
%       10        5        7.83e-8           5   9.58e-5
%       11        4        4.54e-7           4   1.69e-4
%       12        1        6.95e-16          1    0         (no subdomains)



% So it is not possible to make them close unless there are not subdomains.

