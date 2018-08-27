% In this example, our subdomains in x-direction are 
% [0, (a+d)*h], [a*h, (a+2d)*h], [a+d, 1]

RoomData
a=9; d=2;                                   % a and a+d are the interfaces                                         % adding one more sub-domain
maxiter=70;                                 % number of iterations
f1=f(:,1:a+d-1); f2=f(:,a+1:a+2*d-1);           % subdomain source terms
f3=f(:,a+d+1:end);
u1=[gg zeros(J,a+d)];                       % zero initial guess
u2=[zeros(J,2*d+1)];
u3=[zeros(J,J-(a+d)+1) gd];                 % but include boundary values
x1=(0:h:(a+d)*h); x2=(a*h:h:(a+2*d)*h);     % subdomain meshes in x
x3=((a+d)*h:h:1);
z1=zeros(1,a+d+1);                          % useful for plotting 
z2=zeros(1,2*d+1); 
z3=zeros(1,J-(a+d)+2);
z4=zeros(1,J+2); 
u=Solve2d(f,eta,0,J+1,gg,gd);               % direct global solve
u=[z4;u;z4];
mesh(x,y,u);
xlabel('x');
ylabel('y');
zlabel('solution');                         % show solution from direct solve
pause
err(1)=norm(u,'fro');                       % initial error starting with 0
for i=1:maxiter
  u1n=Solve2d(f1,eta,0,a+d,gg,u2(:,d+1));   % solve left subdomain
  u2n=Solve2d(f2,eta,a,a+2*d,u1(:,end-d),u3(:,d+1)); % solve the domain in the middle
  u3n=Solve2d(f3,eta,a+d, J+1 ,u1(:,end),gd);  % solve right subdomain
  u1=u1n;
  u2=u2n;
  u3=u3n;
  ufin=[u1(:,1:a),(u1(:,a+1:a+d)+u2(:,1:d))/2,(u2(:,d+1:end)+u3(:,1:d+1))/2 ,u3(:,d+2:end)];           
  err(i+1)=norm(u(2:end-1,:)-ufin,'fro')  
  mesh(x1,y,[z1;u1n;z1]); hold on           % plot the two iterates
  mesh(x2,y,[z2;u2n;z2]); hold on
  mesh(x3,y,[z3;u3n;z3]); hold off
  xlabel('x');ylabel('y');zlabel('Schwarz iterates');
  pause
end
figure(5)
semilogy((1:maxiter+1)',err)                % plot error decay
xlabel('Iterations');
ylabel('Error');


%% Observations and Analysis

% As the number of subdomains increase, the number of iterations needed for
% convergence become larger.
