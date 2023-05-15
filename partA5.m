
clear all; 
close all;
clc;

%---A1---%

%Specify the constants of the case
T = 2*(10 ^ (-2));
over = 2*10; 
Ts = T/over;
A = 4;
a = 0.5;
Ns=4096;
Fs=1/Ts;

N=100;
bs = (sign (randn(N,1)) + 1)/2;
x = bits_to_2PAM(bs);

%Generate srrc pulse
F=-Fs/2:Fs/Ns:Fs/2-Fs/Ns;
[ph,t] = srrc_pulse(T,over,A,a);
phF=fftshift(fft(ph,Ns).*Ts);
xN=(1/Ts)*upsample(x,over);
t_xN=(0:Ts:N/(1/T)-Ts);


%---A5---%

%Compute Convolution and create time axis
x_conv=conv(ph,xN).*Ts;
t_conv=t(1)+t_xN(1):Ts:t(end)+t_xN(end);

len_total=length(t_conv)*Ts;

%Fourier Transform and periodgram generation
Fx=fftshift(fft(x_conv,Ns)*Ts);
Px=(abs(Fx).^2)/len_total;

k=500;
%Experiment for k repetitions
for i=1:k
    b = (sign(randn(N,1))+1)/2;
    x_test = bits_to_2PAM(b);
    xn=(1/Ts)*upsample(x_test,over);
    x_conv_test=conv(ph,xn)*Ts;
    Fx_test=fftshift(fft(x_conv_test,Ns)*Ts);
    Px=(abs(Fx_test).^2)/len_total;
    P(i,:)=Px;
end


Sx_est=sum(P)/k;
%Compute theoretical Power Spectral Density
Sx_th=(var(x)).*abs((phF).^2)./T;
%Show theoretical and experimental results on same plot
figure()
semilogy(F,Sx_est,'red',F,Sx_th,'blue');
xlabel('Frequency (Hz)');
title('Theoretical and Estimated PSD with 2-PAM configuration for T=2T');
grid on;
legend('S_{X,est} ', 'S_{X,th}')

