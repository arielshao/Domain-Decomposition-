RoomData_Aili                               % add a back door
                                            % larger heater
a=9; d=2;                                   % a and a+d are the interfaces  
maxiter=20;                                 % number of iterations
f1=f(:,1:a+d-1); f2=f(:,a+1:end);           % subdomain source terms
u1=[gg zeros(J,a+d)];                       % zero initial guess
u2=[zeros(J,J-a+1) gd];                     % but include boundary values
x1=(0:h:(a+d)*h); x2=(a*h:h:1);             % subdomain meshes in x 
z1=zeros(1,a+d+1);                          % useful for plotting 
z2=zeros(1,J-a+2); 
z3=zeros(1,J+2); 
u=Solve2d(f,eta,0,J+1,gg,gd);               % direct global solve
u=[z3;u;z3];
mesh(x,y,u);
xlabel('x');
ylabel('y');
zlabel('solution');                         % show solution from direct solve
%pause
err(1)=norm(u,'fro');% return the Frobenius norm                       
% initial error starting with 0
for i=1:maxiter
  u1n=Solve2d(f1,eta,0,a+d,gg,u2(:,d+1));   % solve left subdomain
  u2n=Solve2d(f2,eta,a,J+1,u1(:,end-d),gd); % solve right subdomain
  u1=u1n;
  u2=u2n;
  ufin=[u1(:,1:a),(u1(:,a+1:a+d+1)+u2(:,1:d+1))/2,u2(:,d+2:end)];           
  err(i+1)=norm(u(2:end-1,:)-ufin,'fro')                                     
  mesh(x1,y,[z1;u1n;z1]); hold on           % plot the two iterates
  mesh(x2,y,[z2;u2n;z2]); hold off
  xlabel('x');ylabel('y');zlabel('Schwarz iterates');
  % pause
end
figure(4)
semilogy((1:maxiter+1)',err)                % plot error decay
xlabel('Iterations');
ylabel('Error');
