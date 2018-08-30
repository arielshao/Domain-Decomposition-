RoomData                                    % include problem parameters
a=9; %d=11;                                  % a and a+d are the interfaces
d=11
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

