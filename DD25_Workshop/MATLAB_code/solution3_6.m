%This example solves on 3 subdomains [0, a], [a,3a],[3a,1]
% First solve a Dirichlet problem on the 2nd subdomain
% Next solve a mixed Dirchlet-Neumann problem on 1st and 3rd subdomains
% Then update the interface data and repeat the procedure

RoomData                                    % include problem parameters
f=[zeros(J,1) f zeros(J,1)];                % Robin solvers compute bc
xi=[0 xi 1];                                % thus need to increase size
pe=1e12;                                    % to emulate Dirichlet condition
a=8;                                        % interface position
maxiter=100;
f1=f(:,1:a+1); f2=f(:,a+2:2*a);   
f3=f(:,2*a+1:end);                         % subdomain source terms
g1=zeros(J,1);
g2=zeros(J,1);
x1=(0:h:a*h); x2=(a*h:h:2*a*h); 
x3=(2*a*h:h:1);                            % subdomain mesh in x
z1=zeros(1,a+1); z2=zeros(1,a+1); 
z3=zeros(1,J-2*a+2);                        % for plotting
u=Solve2dR(f,eta,0,J+1,gg*pe,pe*gd,pe,pe);  % global solve
err(1)=norm(u,'fro');                       % initial error starting with 0
th=0.7;                                     % relaxation parameter
e=ones(J,1);                                % construct normal derivative
Na=[speye(J) -spdiags([-e (eta*h^2+4)*e -e]/2,[-1 0 1],J,J)]/h;
Nb=[-spdiags([-e (eta*h^2+4)*e -e]/2,[-1 0 1],J,J) speye(J)]/h;
u1=zeros(J,a+1);
u2=zeros(J,a+1);
u3=zeros(J,J-2*a+1);
figure(3);clf;
for i=1:maxiter
  u2=Solve2d(f2,eta,a,2*a,g1,g2);            % solve dirichlet problem on the 2nd domain
  tb=Nb*[u2(:,1);u2(:,2)]+f1(:,end)*h/2;
  ta=Na*[u2(:,end-1);u2(:,end)]+f3(:,1)*h/2; % extract interface data
  u1=Solve2dR(f1,eta,0,a,gg*pe,tb,pe, 0);           
  u3=Solve2dR(f3,eta,2*a,J+1,ta,gd*pe,0,pe);    % solve dirichlet-neumann problem on 1st and 3rd domain
  g1=th*u1(:,end)+(1-th)*g1;                 % update the interface data
  g2=th*u3(:,1)+(1-th)*g2; 
  ufin=[u1(:,1:a),(u1(:,a+1)+u2(:,1))/2,u2(:,2:end-1),(u2(:,end)+u3(:,1))/2, u3(:,2:end)];
  err(i+1)=norm(u-ufin,'fro');
  mesh(x1,y,[z1;u1;z1]); hold on;           % plot the two iterates
  mesh(x2,y,[z2;u2;z2]); hold on;
  mesh(x3,y,[z3;u3;z3]); hold off;
  xlabel('x');ylabel('y');zlabel('Dirichlet-Neumann iterates');
  %pause
end
figure(4);clf;
semilogy(1:maxiter+1,err)                   % plot error decay
xlabel('Iterations');
ylabel('Error');

%% Observations and Analysis

% a=8  th=0.5 

%                     2  subdomains      vs             3 subdomains
%-----------------------------------------------------------------------
% # of iterations          14                               32
% error                    3.09e-11                         5e-15


% As the number of subdomains increases, more iterations are needed for
% convergence, but our approximate solution converges to a more accurate
% solution.
    