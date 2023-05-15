clear all;
close all;
clc;

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


%---A4---%
K = 500;
N = 100;

%Create n bits series
b = (sign(randn(N,1))+1)/2;

%Encode the bits to the required standard
X = bits_to_4PAM(b);
X_delta=1/Ts*upsample(X,over);
T_plot = 0:Ts:N/2-Ts;
X_delta_conv =conv(X_delta,fi)*Ts;
t_conv = linspace(T_plot(1)+t(1), T_plot(end)+t(end),length(X_delta_conv));

t_all = length(t_conv)*T;

%Calculate the fourier transform of X
XF=fftshift(fft(X_delta_conv,Nf))*Ts; 

%Periodgram of single implementation of X
PxF= (abs(XF).^2)/t_all;

%Sx theoretical
Sx_th = (var(X_delta_conv)/T).*(spect);

%Periodgrams of multiple implementations of X
PxF_lots_of = zeros(K,Nf);

for k = 1:K
    b=(sign(randn(N,1))+1)/2;
    X4 = bits_to_4PAM(b);
    X4_delta =1/Ts*upsample(X4,over);
    X4_delta_conv = conv(X4_delta,fi)*Ts;
    XF4=fftshift(fft(X4_delta_conv,Nf));
    PxF_lots_of(k,:)=(abs(XF4).^2)./t_all;
end
  
%Sx estimated 
Sx_est = sum(PxF_lots_of,1)*Ts./K; 


%semilogy of Sx_est and Sx_th
figure(6)
semilogy(F,Sx_est,'r',F,Sx_th,'b')
title('Theoretical and Estimated PSD with 4-PAM configuration');
xlabel('F');
ylabel('S_X(F)');
legend('S_{X,est} ', 'S_{X,th}')
grid on;
