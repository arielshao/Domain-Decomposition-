a=0;b=1;                                       % domain size
gg=0;gd=1;                                     % boundary conditions
eta=0;                                         % problem parameter (eta-Delta)u=f
%I=2^3;
I=2^4;                                        % number of subdomains, powers of 2
J=2^7-1;                                       % interior mesh points
f=zeros(J,1);                                  % source term zero
Ii=(J+1)/I*(0:I);                              % non-overlapping interface location 
d=4;                                           % 2d is the overlap
maxiter=10;
% figure(5);clf;
figure(9);clf;
[u,err]=SolveDD(f,eta,a,b,gg,gd,Ii,d,zeros(J,1),0,maxiter);   % without coarse
%figure(6);clf;
figure(10);clf;
[u,errc]=SolveDD(f,eta,a,b,gg,gd,Ii,d,zeros(J,1),1,maxiter);  % with MG coarse 
%figure(7);clf;
figure(11);clf;
[u,errcopt]=SolveDD(f,eta,a,b,gg,gd,Ii,d,zeros(J,1),2,2);     % with opt coarse

% figure(8);clf;
figure(12);clf;
semilogy(0:maxiter,err,'-o',0:maxiter,errc,'-+',0:2,errcopt,'-*');
legend('without coarse','with MG coarse','with opt coarse','Location','southeast')
xlabel('Iterations');
ylabel('Error');


%% Observations and Analysis 

% For the original parallel Schwarz iterations, as the number of subdomains
% increases, the rate of convergence decreases.

% For the coarse grid case, as the number of subdomains doubles each time, the 
% rate of convergence also doubles.

% For the optimal coarse grid case, the method always converges to the
% exact solution within two iterations.


% Thus, the parallel Schwarz method with a multigrid type coarse grid is
% scalable.


