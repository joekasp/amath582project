x=csvread('h2o_x_RealTime_Dipole.csv',1,0);
y=csvread('h2o_y_RealTime_Dipole.csv',1,0);
z=csvread('h2o_z_RealTime_Dipole.csv',1,0);

 

 
nn=2000;
t=x(1:nn,1);
m=fix(nn/5);          %simulation period

T=t(nn);        %simulation time
dt=t(2)-t(1);
w=(pi/T)*[0:nn-1];
kick=0.0001;

 

 
pxt=x(1:m,3);
pyt=y(1:m,4);
pzt=z(1:m,5);

 
damp_const=150;
damp = exp(-(t(1:m)-t(1))/damp_const);

 
px=pxt-pxt(1);
py=pyt-pyt(1);
pz=pzt-pzt(1);
px=px.*damp;
py=py.*damp;
pz=pz.*damp;

 
%change to atomic unit
px=0.393456*px;
py=0.393456*py;
pz=0.393456*py;

 

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

 
sigma=(4*pi/(3*137*kick))*w'.*(ax+ay+az);

 
w=27.2114*w;
number=find(w<40);
np=length(number);
figure (2)
plot(w(1:np),sigma(1:np))