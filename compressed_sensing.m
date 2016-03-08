%Compressed Sensing
%AMATH 582 Final Project
%J. Kasper, C. Houferak, J. Radler, S. Sun

%% Step 1: Known Signal
t = linspace(0,1,1000);
S = sin(3*t)+sin(10*t)+sin(40*t); %idealized signal

plot(t,S)
xlabel('Time')
ylabel('Signal Amplitude')
title('Idealized Signal')

n=50; %number of points to sample
samp = S(1:n);

deltaOmega = 1; %frequency difference

%generate Fourier matrix
for j = 1:3
    for k = 1:1000
        F(j,k) = (2/pi)*deltaOmega*sin(t(j)*t(k));
    end
end

cvx_begin;
    variable rec(t);
    minimize(norm(rec,1));
    subject to
    F*rec == samp;

cvx_end;

%% Step 2:
