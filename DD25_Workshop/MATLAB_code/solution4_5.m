RoomData3                                % include problem parameters
pe=1e12;                                    % penalty to emulate Dirichlet
a=9;                                        % interface
maxiter=100;                                 
f1=f(:,1:a-1); f2=f(:,a+1:2*a-1); 
f3=f(:,2*a+1:end);                          % subdomain source terms
x1=(0:h:a*h); x2=(a*h:h:2*a*h);
x3=(2*a*h:h:1);                             % subdomain meshes in x
fn=zeros(J,J+2);                            % Neumann source terms
f1n=fn(:,1:a+1); f2n=fn(:,a+1:2*a+1);
f3n=fn(:,2*a+1:end);
g1=zeros(J,1);
g2=zeros(J,1);                              % updated trace
gzero=zeros(J,1);                           % bc for phi1 and phi2
z1=zeros(1,a+1); 
z2=zeros(1,a+1);                            % for plotting
z3=zeros(1,J-2*a+2);          
z4=zeros(1,J+2);                           
th=0.25;                                     % relaxation parameter
u=Solve2d(f,eta,0,J+1,gg,gd);               % global solve
u=[z4;u;z4];
err(1)=norm(u,'fro');
e=ones(J,1);                                % construct normal derivatives
Na=[speye(J) -spdiags([-e (eta*h^2+4)*e -e]/2,[-1 0 1],J,J)]/h;
Nb=[-spdiags([-e (eta*h^2+4)*e -e]/2,[-1 0 1],J,J) speye(J)]/h;

 for i=1:maxiter
   u1=Solve2d(f1,eta,0,a,gg,g1);              % Dirichlet solves             
   u2=Solve2d(f2,eta,a,2*a,g1,g2);
   u3=Solve2d(f3,eta,2*a,J+1,g2,gd);
   gn1=-Na*[u1(:,end-1);u1(:,end)]-Nb*[u2(:,1);u2(:,2)]-f(:,a)*h/2-f(:,a)*h/2;
   gn2=-Na*[u2(:,end-1);u2(:,end)]-Nb*[u3(:,1);u3(:,2)]-f(:,2*a)*h/2-f(:,2*a)*h/2;
   phi1=Solve2dR(f1n,eta,0,a,pe*gzero,gn1,pe,0); % Neumann solves
   phi2=Solve2dR(f2n,eta,a,2*a,gn1,gn2,0,0);
   phi3=Solve2dR(f3n,eta,2*a,J+1,gn2,pe*gzero,0,pe);
   g1=g1-th*(phi1(:,end)+phi2(:,1));           % update trace
   g2=g2-th*(phi2(:,end)+phi3(:,1));
   ufin=[u1(:,1:end-1),(u1(:,end)+u2(:,1))/2,u2(:,2:end-1),(u2(:,end)+u3(:,1))/2,u3(:,2:end)]; 
   err(i+1)=norm(u(2:end-1,:)-ufin,'fro');
   mesh(x1,y,[z1;u1;z1]); hold on;           % plot Dirichlet iterates
   mesh(x2,y,[z2;u2;z2]); hold on;
   mesh(x3,y,[z3;u3;z3]); hold off;
   xlabel('x');ylabel('y');zlabel('Neumann-Neumann iterates');
% %   pause
 end
semilogy(1:maxiter+1,err)                % plot error decay
xlabel('Iterations');
ylabel('Error');
%% Observations and Analysis 
% a=8

%       th            # of iterations        error 
%   ---------------------------------------------------
%      0.1                64               1.20e-14
%      0.2                23               1.30e-14
%      0.3                50               1.30e-14  
%      0.4                not converge    
%      0.5                diverge


% As the number of subdomains increases, the number of iterations needed
% for convergence changes; it takes less iterations to converge to a more
% accuarate solution for ttheta=0.1 and theta=0.2, but it no longer
% converges for theta=0.4.
