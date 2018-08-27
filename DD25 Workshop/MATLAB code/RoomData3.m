eta=0; J=19;                                 % no of interior mesh points
h=1/(J+1);                                   % mesh size 
x=(0:h:1); y=x;                              % mesh in x and y direction
f=zeros(J,J);                                % interior source function
xi=x(2:end-1); yi=xi;                        % interior mesh in x-y direction 
f([yi>0.4 & yi<0.6],[xi>0.4 & xi<0.6])=50;   % heater placement
% f([yi>0.1 & yi<0.9],[xi>0.1 & xi<0.9])=50;   % floor heater
gg=0.3*ones(J,1); gg(yi>0.5 & yi<0.9)=1;     % door right hand side
gd=zeros(J,1);