x=csvread('h2o_x_RealTime_Dipole.csv',1,0);
y=csvread('h2o_y_RealTime_Dipole.csv',1,0);
z=csvread('h2o_z_RealTime_Dipole.csv',1,0);


nn=length(x(:,1));
t=x(:,1);
m=50;          %simulation period
n=250;         %ideal period  
T=t(end);        %simulation time
dt=t(2)-t(1);
w=(2*pi/T)*[0:nn-1];
kick=0.0001;


pxt=x(:,3);
pyt=y(:,4);
pzt=z(:,5);

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

%%

FF=idst(eye(n,n));
F=FF(1:m,:);
cvx_begin;
    variable ax(n,1);
    minimize( norm(ax,1) );
    subject to 
    F*ax==px;
cvx_end;

cvx_begin;
    variable ay(n,1);
    minimize( norm(ay,1) );
    subject to 
    F*ay==py;
cvx_end;

cvx_begin;
    variable az(n,1);
    minimize( norm(az,1) );
    subject to 
    F*az==pz;
cvx_end;


s=(4*pi/(3*kick*137))*w'.*(ax+ay+az);

%%
figure (1)
plot(w(1:20),s(1:20))

%%
% pxt=x(1:n,3);
% pyt=y(1:n,4);
% pzt=z(1:n,5);
% 
% px=pxt-pxt(1);
% py=pyt-pyt(1);
% pz=pzt-pzt(1);

ax=fftshift(fft(px));
ay=fftshift(fft(py));
az=fftshift(fft(pz));

sigma=(4*pi/(3*137*kick))*w'.*imag(ax+ay+az);


%%

% sigma=fftshift(sigma);
figure (2)
plot(w(1:nn/2)*27.2114,sigma(1:nn/2))

