a=0;b=1;                                       % domain size
gg=0;gd=1;                                     % boundary conditions
%eta=0;                                         % problem parameter (eta-Delta)u=f
eta=2;
I=2^3;                                         % number of subdomains, powers of 2
J=2^7-1;                                       % interior mesh points
f=zeros(J,1);                                  % source term zero
Ii=(J+1)/I*(0:I);                              % non-overlapping interface location 
d=4;                                           % 2d is the overlap
maxiter=10;
[u,err]=SolveDD(f,eta,a,b,gg,gd,Ii,d,zeros(J,1),0,maxiter);   % without coarse
[u,errc]=SolveDD(f,eta,a,b,gg,gd,Ii,d,zeros(J,1),1,maxiter);  % with MG coarse 
[u,errcopt]=SolveDD(f,eta,a,b,gg,gd,Ii,d,zeros(J,1),2,2);     % with opt coarse

semilogy(0:maxiter,err,'-o',0:maxiter,errc,'-+',0:2,errcopt,'-*');
legend('without coarse','with MG coarse','with opt coarse','Location','southeast')
xlabel('Iterations');
ylabel('Error');

%% Observations and Analysis
% If we take eta=5 instead of 0, the present optimal coarse space in the
% code is not optimal anymore.

% The choice of the optimal coarse space depends on the elliptic
% operator. The optimal coarse space is this code is designed for the
% Poisson's equation only.
