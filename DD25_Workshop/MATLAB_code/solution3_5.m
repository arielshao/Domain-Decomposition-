%This code simulates an insulated wall on the right edge of the domain;
% that is, partial u/ partial n=0 for x=1;

RoomData                                    % include problem parameters
f=[zeros(J,1) f zeros(J,1)];                % Robin solvers compute bc
xi=[0 xi 1];                                % thus need to increase size
pe=1e12;                                    % to emulate Dirichlet condition
a=8;                                        % interface position
maxiter=200;
f1=f(:,2:a); f2=f(:,a+1:end);               % subdomain source terms
g=zeros(J,1);
x1=(0:h:a*h); x2=(a*h:h:1);                 % subdomain mesh in x
z1=zeros(1,a+1); z2=zeros(1,J-a+2);         % for plotting
u=Solve2dR(f,eta,0,J+1,gg*pe,gd,pe,0);      % global solve
err(1)=norm(u,'fro');                       % initial error starting with 0
th=0.7;                                     % relaxation parameter
e=ones(J,1);                                % construct normal derivative
Na=[speye(J) -spdiags([-e (eta*h^2+4)*e -e]/2,[-1 0 1],J,J)]/h;
figure(100);clf;
for i=1:maxiter
  u1=Solve2d(f1,eta,0,a,gg,g);              % solve left subdomain
  ta=Na*[u1(:,end-1);u1(:,end)]+f2(:,1)*h/2;% extract interface data
  u2=Solve2dR(f2,eta,a,J+1,ta,gd,0,0);     % solve right subdomain
  g=th*u2(:,1)+(1-th)*g;                    % relax Dirichlet trace 
  ufin=[u1(:,1:a),(u1(:,a+1)+u2(:,1))/2,u2(:,2:end)];
  err(i+1)=norm(u-ufin,'fro');
  mesh(x1,y,[z1;u1;z1]); hold on;           % plot the two iterates
  mesh(x2,y,[z2;u2;z2]); hold off;
  xlabel('x');ylabel('y');zlabel('Dirichlet-Neumann iterates');
  pause
end
figure(101);clf;
semilogy(1:maxiter+1,err)                   % plot error decay
xlabel('Iterations');
ylabel('Error');


%% Observations and Analysis 

% At iteration 20:

%           Homogeneous Neumann      vs    Homegeneous Dirichlet 
%-------------------------------------------------------------------
%  error       1.52e-4                         1.06e-5
% solution  monotonically decreasing         
%          to zero on the right domain         like a saddle
%-------------------------------------------------------------------
% The insulated wall case is slow in convergence.
% Graphically, it gives a more resonable simulation of the heating of a
% room.

