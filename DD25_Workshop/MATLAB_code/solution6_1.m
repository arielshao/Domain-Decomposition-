
a=0;b=1;                                       % domain size
gg=0;gd=1;                                     % boundary conditions
eta=0;                                         % problem parameter (eta-Delta)u=f
I=2^2;                                         % number of subdomains, powers of 2
J=2^7-1;                                       % interior mesh points
f=zeros(J,1);                                  % source term zero
Ii=(J+1)/I*(0:I);                              % non-overlapping interface location 
d=4;                                           % 2d is the overlap
maxiter=10;
figure(1);clf;
[u,err]=SolveDD(f,eta,a,b,gg,gd,Ii,d,zeros(J,1),0,maxiter);   % without coarse
figure(2);clf;
[u,errc]=SolveDD(f,eta,a,b,gg,gd,Ii,d,zeros(J,1),1,maxiter);  % with MG coarse 
figure(3);clf;
[u,errcopt]=SolveDD(f,eta,a,b,gg,gd,Ii,d,zeros(J,1),2,2);     % with opt coarse

figure(4);clf;
semilogy(0:maxiter,err,'-o',0:maxiter,errc,'-+',0:2,errcopt,'-*');
legend('without coarse','with MG coarse','with opt coarse','Location','southeast')
xlabel('Iterations');
ylabel('Error');

%% Observations and Analysis 

% The original parallel Schwarz method is slow in convergence; it has not
% converged to the true solution after 10 iterations.

% The parallel Schwarz with a multigrid type coarse grid is slightly faster
% in convergence; it converges to the exact solution with an error of 9.05e-4 after
% 10 iterations

% The parallel Schwarz with an optimal coarse grid is extremely fast in
% convergence; it converges to the exact solution with an error of 6e-15 within 2
% iterations.