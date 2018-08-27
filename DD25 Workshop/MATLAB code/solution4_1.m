RoomData                                    % include problem parameters
pe=1e12;                                    % penalty to emulate Dirichlet
a=8;                                        % interface
maxiter=20;                                 
f1=f(:,1:a-1); f2=f(:,a+1:end);             % subdomain source terms
x1=(0:h:a*h); x2=(a*h:h:1);                 % subdomain meshes in x
fn=zeros(J,J+2);                            % Neumann source terms
f1n=fn(:,1:a+1); f2n=fn(:,a+1:end);
g=zeros(J,1);                               % updated trace
gzero=zeros(J,1);                           % bc for phi1 and phi2
z1=zeros(1,a+1);                            % for plotting
z2=zeros(1,J-a+2);          
z3=zeros(1,J+2);                            
th=0.4;                                     % relaxation parameter
u=Solve2d(f,eta,0,J+1,gg,gd);               % global solve
u=[z3;u;z3];
err(1)=norm(u,'fro');
e=ones(J,1);                                % construct normal derivatives
Na=[speye(J) -spdiags([-e (eta*h^2+4)*e -e]/2,[-1 0 1],J,J)]/h;
Nb=[-spdiags([-e (eta*h^2+4)*e -e]/2,[-1 0 1],J,J) speye(J)]/h;

for i=1:maxiter
  u1=Solve2d(f1,eta,0,a,gg,g);              % Dirichlet solves             
  u2=Solve2d(f2,eta,a,J+1,g,gd);
  gn=-Na*[u1(:,end-1);u1(:,end)]-Nb*[u2(:,1);u2(:,2)]-f(:,a)*h/2-f(:,a)*h/2;
  phi1=Solve2dR(f1n,eta,0,a,pe*gzero,gn,pe,0); % Neumann solves
  phi2=Solve2dR(f2n,eta,a,J+1,gn,pe*gzero,0,pe);
  g=g-th*(phi1(:,end)+phi2(:,1));           % update trace
  ufin=[u1(:,1:end-1),(u1(:,end)+u2(:,1))/2,u2(:,2:end)];
  err(i+1)=norm(u(2:end-1,:)-ufin,'fro');
  mesh(x1,y,[z1;u1;z1]); hold on;           % plot Dirichlet iterates
  mesh(x2,y,[z2;u2;z2]); hold off;
  xlabel('x');ylabel('y');zlabel('Neumann-Neumann iterates');
  pause
end
semilogy(1:maxiter+1,err)                % plot error decay
xlabel('Iterations');
ylabel('Error');

%% Obeservations and Analysis 
% As the number of iterations increases, the Neumann-Neumann method
% converges to the global solution.

% The iterations alternate between super and sub solutions of the global
% solution; that is, one approximation is greater than the true solution
% and then next approximation is smaller than the true solution.



