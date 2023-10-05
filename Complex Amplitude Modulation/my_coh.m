function z = my_coh(y,fc)
% z = extracted demodulated signal
% y = modulated signal
% fc = carrier frequency
Ac = input('enter carrier amplitude: ');
phi = input('enter phase difference of local oscillator: ');
fs = input('enter sampling frequency: ');
ts = input('enter time length: ');
fcut = fc/(fs/2);
tiv =1/fs;
t = 0:tiv:ts-tiv;
lo = Ac*cos(2*pi*fc*t+phi);
v = y.*lo; % product modulator
% passing into my_lpf
z = my_lpf(v, fs, fcut);

figure;
subplot(2,1,1), plot(t, y), title('modulated signal');
subplot(2,1,2), plot(t, z), title('recovered signal');

end
