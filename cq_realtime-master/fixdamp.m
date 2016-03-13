sample=5:5:50000; 
sample_cut=(sample(1:500))';
x=csvread('h2o_x_RealTime_Dipole.csv',1,0);
y=csvread('h2o_y_RealTime_Dipole.csv',1,0);
z=csvread('h2o_z_RealTime_Dipole.csv',1,0);

 



t=x(sample(1:2000),1);
m=length(sample)/5;          %simulation period
nn=m*5;

T=t(end);        %simulation time
dt=t(2)-t(1);
w=(pi/T)*[1:m];
kick=0.0001;

 

 
pxt=x(1:m,3);
pyt=y(1:m,4);
pzt=z(1:m,5);

 
damp_const=200;
damp = exp(-(t(1:m)-t(1))/damp_const);

 
px=pxt-pxt(1);
py=pyt-pyt(1);
pz=pzt-pzt(1);
px=px.*damp;
py=py.*damp;
pz=pz.*damp;

%% 
%change to atomic unit
px=0.393456*px;
py=0.393456*py;
pz=0.393456*py;

 
%%
FF=idst(eye(nn,nn));
F=FF(1:m,:);
 
cvx_begin;
    variable ax(nn,1);
    minimize(norm(ax,1));
    subject to
    norm((F*ax-px))<0.00001;
cvx_end;

cvx_begin;
    variable ay(nn,1);
    minimize(norm(ay,1));
    subject to
    norm((F*ay-py))<0.00001;
cvx_end;

cvx_begin;
    variable az(nn,1);
    minimize(norm(az,1));
    subject to
    norm((F*az-pz))<0.00001;
cvx_end;

%%

ax=fft(px);
ay=fft(py);
az=fft(pz);

sigma=(4*pi/(3*137*kick))*w'.*(ax+ay+az);

 
w=27.2114*w;
number=find(w<40);
np=length(number);
figure (2)
plot(w(1:np),sigma(1:np))