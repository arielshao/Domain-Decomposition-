RoomData                                    % include problem parameters
a=10; %d=11;                                  % a and a+d are the interfaces
d=0;
pe=1e12;

maxiter=20;

f=[zeros(J,1) f zeros(J,1)]; 
f1=f(:,2:a); f2=f(:,a+1:end);           % subdomain source terms
z=zeros(J,1);

G1=zeros(a+d+1,1);G1(end-d)=1;              % trace operators
G2=zeros(J-a+2,1);G2(d+1)=1;

b=[Solve2d(f1,eta,0,a+d,gg,z)*G1;           % construct interface system    
   Solve2dR(f2,eta,a,J+1,z,gd,pe,pe)*G2];          % T*g=b, g unknowns at interfaces

g=zeros(2*J,1); 
th=0.5;                                     % relaxation parameter
e=ones(J,1);                                % construct normal derivative
Na=[speye(J) -spdiags([-e (eta*h^2+4)*e -e]/2,[-1 0 1],J,J)]/h;

u1=Solve2d(f1,eta,0,a,gg,g(J+1:2*J));  
ta=Na*[u1(:,end-1);u1(:,end)]+f2(:,1)*h/2;

T=@(g) [g(1:J)-Solve2d(zeros(size(f1,1),size(f1,2)),eta,0,a+d,gg,g(J+1:2*J))*G1;
        g(J+1:2*J)-Solve2dR(zeros(size(f2,1),size(f2,2)),eta,a,J+1,ta,gd,0,pe)*G2];
                            % zero initial guess
ri(1)=norm(b-T(g));                         % initial residual
for i=1:maxiter 
   g(J+1:2*J)=th*g(1:J)+(1-th)*g(J+1:2*J);
   g=g-T(g)+b; 
   ri(i+1)=norm(b-T(g));                     % keep residual for plotting
end
[g,fl,r,it,rk]=gmres(T,b);                  % use GMRES Krylov solver
semilogy(0:length(ri)-1,ri,'-o',0:length(rk)-1,rk,'-+');

xlabel('iteration')
ylabel('residual')
legend('Iterative','Krylov')

