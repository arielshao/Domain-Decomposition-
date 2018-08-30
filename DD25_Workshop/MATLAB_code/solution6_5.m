% a=0;b=1;                                       % domain size
a=0;b=2;                                   % increase the domain size
ddx=0.25;                                  % each subdomain size
gg=0;gd=1;                                 % boundary conditions
eta=0;                                     % problem parameter (eta-Delta)u=f
% I=2^2;               % number of subdomains proportional to domain size
I=(b-a)/ddx;
J=(2^7)*(b-a)-1;   % increase interior mesh points proportionally
f=zeros(J,1);                                  % source term zero
Ii=(J+1)/I*(0:I);                              % non-overlapping interface location 
d=4;                                           % 2d is the overlap
maxiter=10;
figure(20);clf;
[u,err]=SolveDD(f,eta,a,b,gg,gd,Ii,d,zeros(J,1),0,maxiter);   % without coarse
[u,errc]=SolveDD(f,eta,a,b,gg,gd,Ii,d,zeros(J,1),1,maxiter);  % with MG coarse 
[u,errcopt]=SolveDD(f,eta,a,b,gg,gd,Ii,d,zeros(J,1),2,2);     % with opt coarse

semilogy(0:maxiter,err,'-o',0:maxiter,errc,'-+',0:2,errcopt,'-*');
legend('without coarse','with MG coarse','with opt coarse','Location','southeast')
xlabel('Iterations');
ylabel('Error');


%% Observations and Analysis 

% Case I: eta=0;
% If we let the domain (a,b) grow proportionally when we add more
% subdomains, the one level method is scalable.


% Case II: eta=5;
% For this case, the one level method is scalable as well.
