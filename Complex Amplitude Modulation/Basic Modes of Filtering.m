% provide n,fs,fcut,ftype
% fs must be atleast 2*fc as per nyquist criteria
fsx = 10e3;
fxcut1 = 60/(fsx/2) ; %frequency normalization to lie between 0 and 1 using Nyquist Criteria
fxcut2 = 200/(fsx/2) ;
tix = 1/fsx;
t = 0:tix:1-tix;
x = square(2*pi*50*t);
%% -------------------------------------------------------------
% Part a
ylo1 = my_lpf(x, fsx, fxcut1);
ylo2 = my_lpf(x, fsx, fxcut2);
figure
subplot(3,1,1), plot(t,x,'r'), xlim([0.4 0.45]), title('input signal');
subplot(3,1,2), plot(t,ylo1,'g'), xlim([0.4 0.45]), title('low filtered signal with fc = 60');
subplot(3,1,3), plot(t,ylo2,'b'), xlim([0.4 0.45]), title('low filtered signal with fc = 200');
figure, freqz(x);
figure, freqz(ylo1);
figure, freqz(ylo2);
%% ----------------------------------------------------------------
% Part b
fxcut = 100/(fsx/2);
yhi = my_hpf(x,fsx, fxcut);
figure, freqz(x);
figure, freqz(yhi);
%% ----------------------------------------------------------------
% Part c
fclo = 100/(fsx/2);
fchi = 200/(fsx/2);
ybp = my_bpf(x,fsx, fclo, fchi);
figure, freqz(x);
figure, freqz(ybp);
