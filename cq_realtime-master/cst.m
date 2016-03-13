x=csvread('h2o_x_RealTime_Dipole.csv',1,0);
y=csvread('h2o_y_RealTime_Dipole.csv',1,0);
z=csvread('h2o_z_RealTime_Dipole.csv',1,0);
% read the data of dipole moment
 


nn=2000;            % number of data points (long edge of F matrix)
t=x(1:nn,1);        % time
m=fix(nn/5);          %number of used data points  (short edge of F matrix)

T=t(nn);           %simulated simulation time. actual simulation time will be T/5.
dt=t(2)-t(1);       % time step
w=(pi/T)*[0:nn-1];  % frequency
kick=0.0001;        % amplitude of electric field E 

 

 
pxt=x(1:m,3);       
pyt=y(1:m,4);
pzt=z(1:m,5);
% pick the data of angular momentum
 
damp_const=1500;
damp = exp(-(t(1:m)-t(1))/damp_const);

 
px=pxt-pxt(1);
py=pyt-pyt(1);
pz=pzt-pzt(1);
% subtract mean value

px=px.*damp;
py=py.*damp;
pz=pz.*damp;

% damp

%change to atomic unit
px=0.393456*px;
py=0.393456*py;
pz=0.393456*py;

 %% 
% do L1 optimization
 
FF=idst(eye(nn,nn));
F=FF(1:m,:);
 
cvx_begin;
    variable ax(nn,1);
    minimize(norm(ax,1));
    subject to
    F*ax==px;
cvx_end;

cvx_begin;
    variable ay(nn,1);
    minimize(norm(ay,1));
    subject to
    F*ay==py;
cvx_end;

cvx_begin;
    variable az(nn,1);
    minimize(norm(az,1));
    subject to
    F*az==pz;
cvx_end;

%%
 
sigma=(4*pi/(3*137*kick))*w'.*(ax+ay+az);  % sigma=blabla*(axx+ayy+azz)

 
w=27.2114*w;                               %change frequency w into eq
number=find(w<40);                         % confine the w to the frequency interval of we want
np=length(number);                         
figure (2)
plot(w(1:np),sigma(1:np))