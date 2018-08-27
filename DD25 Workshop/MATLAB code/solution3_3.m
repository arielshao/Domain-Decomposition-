RoomData                                    % include problem parameters
f=[zeros(J,1) f zeros(J,1)];                % Robin solvers compute bc
xi=[0 xi 1];                                % thus need to increase size
pe=1e12;                                    % to emulate Dirichlet condition
a=4;                                        % interface position
maxiter=50;
f1=f(:,2:a); f2=f(:,a+1:end);               % subdomain source terms
g=zeros(J,1);
x1=(0:h:a*h); x2=(a*h:h:1);                 % subdomain mesh in x
z1=zeros(1,a+1); z2=zeros(1,J-a+2);         % for plotting
u=Solve2dR(f,eta,0,J+1,gg*pe,pe*gd,pe,pe);  % global solve
err(1)=norm(u,'fro');                       % initial error starting with 0
th=(1-a*h)*2/(a*h+3*(1-a*h));               % relaxation parameter
% th= (1-a*h)
% th=0.553
e=ones(J,1);                                % construct normal derivative
Na=[speye(J) -spdiags([-e (eta*h^2+4)*e -e]/2,[-1 0 1],J,J)]/h;
for i=1:maxiter
  u1=Solve2d(f1,eta,0,a,gg,g);              % solve left subdomain
  ta=Na*[u1(:,end-1);u1(:,end)]+f2(:,1)*h/2;% extract interface data
  u2=Solve2dR(f2,eta,a,J+1,ta,gd,0,pe);     % solve right subdomain
  g=th*u2(:,1)+(1-th)*g;                    % relax Dirichlet trace 
  ufin=[u1(:,1:a),(u1(:,a+1)+u2(:,1))/2,u2(:,2:end)];
  err(i+1)=norm(u-ufin,'fro');
  mesh(x1,y,[z1;u1;z1]); hold on;           % plot the two iterates
  mesh(x2,y,[z2;u2;z2]); hold off;
  xlabel('x');ylabel('y');zlabel('Dirichlet-Neumann iterates');
 % pause
end
semilogy(1:maxiter+1,err)                   % plot error decay
xlabel('Iterations');
ylabel('Error');

%% Observations and Analysis 
% a=8

%       th            # of iterations        error 
%   ---------------------------------------------------
%      0.1               120               3.09e-11
%      0.2                53               3.09e-11
%      0.3                30               3.09e-11  
%      0.4                18               3.09e-11
%      0.5                14               3.09e-11 
%      0.6                26               3.09e-11 
%      0.7                14               3.09e-11
%      0.8                86               3.09e-11
%      0.9                472              3.09e-11

% Both th=0.5 and th=0.7 are optimal values in this case.
% The optimal value given in the lectures is th_opt=2b/(a+3b)=0.553
% The optimal value th=0.5 found heuristically is quite close to the
% theoretical optimal value theta.

% If we take th=(1-a*h)*2/(a*h+3*(1-a*h))(the optimal value of theta), then
% the method should always converges. But in practice, it still diverges
% for a=2 and a=3.

