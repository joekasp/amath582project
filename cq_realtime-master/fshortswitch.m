x=csvread('h2o_x_RealTime_Dipole.csv',1,0);
y=csvread('h2o_y_RealTime_Dipole.csv',1,0);
z=csvread('h2o_z_RealTime_Dipole.csv',1,0);

 

 
% nn=fix(length(x(:,1))/20);
nn=2000;
t=x(1:nn,1);
        %simulation period
n=length(x(:,1));         %ideal period  
T=t(nn);        %simulation time
dt=t(2)-t(1);
w=(2*pi/T)*[0:nn-1];
kick=0.0001;

 

 
pxt=x(1:nn,3);
pyt=y(1:nn,4);
pzt=z(1:nn,5);

 
damp_const=150;
damp = exp(-(t-t(1))/damp_const);

 
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

 p=px+py+pz;

 
a=fft(p);


sigma=-(4*pi/(3*137*kick))*w'.*imag(a);

 
w=27.2114*w;
number=find(w<40);
np=length(number);
figure (2)
plot(w(1:np),sigma(1:np))