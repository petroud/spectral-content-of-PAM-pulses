clear all;
close all;
clc;

%----------------------------%
% Author: Petrou Dimitrios
% Year: 2023  
% TU of Crete
% Telecommunication Systems I
%----------------------------%


                         %---Part A---%
%----------------------------------------------------------------%
%---A1---%
T = 10.^-2;
over = 10;
Ts = T/over;
A = 4;
a = 0.5;

[fi, t] = srrc_pulse(T, over, A, a);
figure(1);
plot(t, fi, 'DisplayName', 'a=0.5');
title('SRRC pulse for A=4 / a=0.5 / T= 10^{-2}s');
xlabel('Time (s)');
ylabel('SRRC Pulse Amplitude');

%Constants of the case
Fs = 1/Ts;
Nf = 2048;

%Defining the frequencies with equal distances
F = (-Fs/2):(Fs/Nf):(Fs/2 - Fs/Nf);

%Perform Fourier transform using fft()
%and use fftshift to rearrange it
ft = fftshift(fft(fi,Nf)*Ts);

%Calculate the energy density spectrum of the Fourier transform
spect = abs(ft).^2;

figure(2)
semilogy(F, spect);
title('Energy Spectral Density |Ö(F)|^2');
xlabel('Frequency (Hz)');
ylabel('ESD');


%---A2---%
N = 100;
%Create n bits series
b = (sign(randn(N,1))+1)/2;

%Encode the bits to the required standard
Xn = bits_to_2PAM(b);
%Time of 
T_plot = 0:Ts:N-Ts;

%Calculate Xä and the convolution with the pulse
X_delta=1/Ts*upsample(Xn,over);
X_delta_conv =conv(X_delta,fi)*Ts;

%Define the convolute time vector
t_conv = linspace(T_plot(1)+t(1), T_plot(end)+t(end),length(X_delta_conv));

figure(3)
plot(t_conv,X_delta_conv)
title('X(t) from 2-PAM');
xlabel('Time in sec');
ylabel('X(t)');
grid on


%---A3---%
t_all = length(t_conv)*T;

%Calculate the fourier transform of X inside the
%Periodgram of single implementation of X
PxF= ((abs(fftshift(fft(X_delta_conv,Nf))).^2)*Ts)/t_all;

%plot 
figure(4);
subplot(2,1,1)
plot(F,PxF)
title('Plot of single implementation periodogram with 2-PAM configuration');
xlabel('F');
ylabel('P_X(F)')
grid on;

%semilogy
subplot(2,1,2)
semilogy(F, PxF);
title('Semilogy of single implementation periodogram with 2-PAM configuration');
xlabel('F');
ylabel('P_X(F)')
grid on;


K = 500;
%Sx theoretical
Sx_th = (var(X_delta_conv)/T).*(spect);

%Periodgrams of multiple implementations of X
PxF_lots_of = zeros(K,Nf);


for k = 1:K
    b=(sign(randn(N,1))+1)/2;
    XR = bits_to_2PAM(b);
    X_deltaR=1/Ts*upsample(XR,over);
    X_delta_convR= conv(X_deltaR,fi)*Ts;
    PxF=((abs(fftshift(fft(X_delta_convR,Nf)))).^2)./t_all;
    PxF_lots_of(k,:)=PxF;
end
  
%Sx estimated 
Sx_est = sum(PxF_lots_of,1)*Ts./K; 

figure(5);
semilogy(F,Sx_est,'r',F,Sx_th,'b')
title('Theoretical and Estimated PSD with 2-PAM configuration');
xlabel('F');
ylabel('S_X(F)');
legend('S_{X,est} ', 'S_{X,th}')
grid on;