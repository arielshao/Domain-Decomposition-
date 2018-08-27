RoomData                                    % include problem parameters
f=[zeros(J,1) f zeros(J,1)];                % Robin solvers compute bc
xi=[0 xi 1];                                % thus need to increase size
pe=1e12;                                    % to emulate Dirichlet condition
a=9; d=2;                                   % a and a+d are the interfaces
maxiter=100;
f1=f(:,1:a+d+1);f2=f(:,a+1:end);            % subdomain source terms
u1=[gg zeros(J,a+d)];                       % zero initial guess
u2=[zeros(J,J-a+1) gd];                     % but include boundary values
x1=(0:h:(a+d)*h); x2=(a*h:h:1);             % subdomain mesh in x
z1=zeros(1,a+d+1);                          % useful for plotting 
z2=zeros(1,J-a+2);
z3=zeros(1,J+2);                          
u=Solve2dR(f,eta,0,J+1,pe*gg,pe*gd,pe,pe);  % direct global solve
u=[z3;u;z3];
mesh(x,y,u);
xlabel('x');
ylabel('y');
zlabel('solution');                         % show solution from direct solve
%pause
err(1)=norm(u,'fro');
% if d>0
%   p=((pi^2+eta)/(d*h))^(1/3);               % optimized choice with overlap
%                                             % but eta=0??
% else
%   p=sqrt(pi^2/h);                           % optimized choice without overlap 
% end

p=4;
e=ones(J,1);                                % construct normal derivatives
Na=[speye(J) -spdiags([-e (eta*h^2+4)*e -e]/2,[-1 0 1],J,J)]/h;
Nb=[-spdiags([-e (eta*h^2+4)*e -e]/2,[-1 0 1],J,J) speye(J)]/h;
for i=1:maxiter                             % extract interface data
  tb=Nb*[u2(:,d+1);u2(:,d+2)]+f2(:,d+1)*h/2+p*u2(:,d+1);
  ta=Na*[u1(:,end-d-1);u1(:,end-d)]+f1(:,end-d)*h/2+p*u1(:,end-d); 
  u1n=Solve2dR(f1,eta,0,a+d,pe*gg,tb,pe,p); % solve left subdomain
  u2n=Solve2dR(f2,eta,a,J+1,ta,gd,p,pe);    % solve right subdomain
  u1=u1n;
  u2=u2n;
  mesh(x1,y,[z1;u1n;z1]); hold on;mesh(x2,y,[z2;u2n;z2]); hold off
  xlabel('x'); ylabel('y'); zlabel('Optimized Schwarz iterates');
  %pause
  ufin=[u1n(:,1:a),(u1n(:,a+1:a+d+1)+u2n(:,1:d+1))/2,u2n(:,d+2:end)];
  err(i+1)=norm(u(2:end-1,:)-ufin,'fro')
end
figure(10);
semilogy((1:maxiter+1)',err)                % plot error decay
xlabel('Iterations');
ylabel('Error');
%% Observations and Analysis 

%       p                     # of iterations        error 
%   ---------------------------------------------------
%      100                          94                3.6e-14
%      50                           80                4.2e-14
%      25                           64                4.4e-14  
%      12                           42                4.6e-14
%       6                           24                4.8e-14 
%       5                           19                4.7e-14 
%       4                           19                4.7e-14
% p=((pi^2+eta)/(d*h))^(1/3)=4.697  17                4.7e-14

% The convergence speed depend on the parameter p.
% The formula in the code is close to the best I can find.
