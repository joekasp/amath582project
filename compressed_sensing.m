%Compressed Sensing
%AMATH 582 Final Project
%J. Kasper, C. Houferak, J. Radler, S. Sun

%% Step 1: Known Signal
clc;
clear;
N = 1000; %total time length
t2 = linspace(0,1,N+1); t=t2(1:N);
S = 4*sin(15*t*(2*pi))+2*sin(7*t*(2*pi))+3*sin(3*t*(2*pi)); %idealized signal

n=200; %number of points to sample
samp = S(1:n)';

figure(1)
subplot(2,1,1)
plot(t,S)
xlabel('Time')
ylabel('Signal Amplitude')
title('Idealized Signal')

subplot(2,1,2)
plot(t(1:n),samp)
xlabel('Time')
ylabel('Signal Amplitude')
title('Short Sample')

deltaOmega = 1; %frequency difference

%generate Fourier matrix
for j = 1:n
    for k = 1:N
        F(j,k) = (2/pi)*deltaOmega*sin(j*t(k));
    end
end

%solve the equation F*x = samp
sigma = 0.05;
opts = spgSetParms('verbosity',0);         % Turn off the SPGL1 log output
x = spg_bpdn(F, samp, sigma, opts);

%norm(samp-F*x)

csSignal = zeros(size(t));
for i = 1:N
    csSignal = csSignal + sin(deltaOmega*t*i)*x(i);
end

x = x/max(x); %normalize

figure(2)
plot(t,csSignal)
xlabel('Time')
ylabel('Signal Amplitude')
title('Compressed Sensing Signal')


St2 = abs(fftshift(fft(S)));
St = St2(N/2+1:N)/max(St2);

sampt2 = abs(fftshift(fft(samp)));
sampt = sampt2(n/2:n)/max(sampt2);

for i = 1:length(t)
    freq(i) = i/(2*pi);
end

for i = 1:length(St)
    freq2(i) = i;
end

for i = 1:length(sampt)
    freq3(i) = (N/n)*(i-1);
end

figure(3)
hold on
subplot(3,1,1), plot(freq2,St);
title('Ideal Frequency Spectrum')
xlabel('Frequency (Hz)')
ylabel('Normalized Intensity')
axis([0 100 -.2 1])
subplot(3,1,2), plot(freq,x)
title('Compressed Sensing Frequency Spectrum')
xlabel('Frequency (Hz)')
ylabel('Normalized Intensity')
axis([0 100 -.2 1])
subplot(3,1,3), plot(freq3,sampt)
title('Fourier Transform Frequency Spectrum')
xlabel('Frequency (Hz)')
ylabel('Normalized Intensity')
axis([0 100 -.2 1])
hold off

%% Step 2:
