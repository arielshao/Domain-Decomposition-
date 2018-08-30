a=0;b=1;                                       % domain size
gg=0;gd=1;                                     % boundary conditions
eta=0;                                         % problem parameter (eta-Delta)u=f
I=2^2;                                         % number of subdomains, powers of 2
% J=2^7-1;                                       % interior mesh points
J=2^8-1;
f=zeros(J,1);                                  % source term zero
Ii=(J+1)/I*(0:I);                              % non-overlapping interface location 
d=4;                                           % 2d is the overlap
maxiter=10;
figure(13);clf;
[u,err]=SolveDD(f,eta,a,b,gg,gd,Ii,d,zeros(J,1),0,maxiter);   % without coarse
figure(14);clf;
[u,errc]=SolveDD(f,eta,a,b,gg,gd,Ii,d,zeros(J,1),1,maxiter);  % with MG coarse 
figure(15);clf;
[u,errcopt]=SolveDD(f,eta,a,b,gg,gd,Ii,d,zeros(J,1),2,2);     % with opt coarse

figure(16);clf;
semilogy(0:maxiter,err,'-o',0:maxiter,errc,'-+',0:2,errcopt,'-*');
legend('without coarse','with MG coarse','with opt coarse','Location','southeast')
xlabel('Iterations');
ylabel('Error');


%% Observations and Analysis

% For the original parallel Schwarz method, as the number of mesh points
% increases, the rate of convergence decreases.

% For the coarse grid case, as the number of mesh points doubles (the mesh
% size is halved), the rate of convergence decreases as well.

% For the optimal coarse grid case, the method always converges within two
% iterations.

% The matrix size of the overlap region becomes larger.

% The parallel schwarz with optimal coarse grid is scalable in this case.