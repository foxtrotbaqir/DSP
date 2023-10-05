function y = my_dsbsc(m, fc)
% where y = output
% m = message signal
% fc = carrier frequency
% Double-sideband suppressed-carrier transmission (DSB-SC) is transmission
% in which frequencies produced by amplitude modulation (AM) are symmetrically
% spaced above and below the carrier frequency and the carrier level is reduced
% to the lowest practical level, ideally being completely suppressed.
Ac = input('Enter the value of carrier signal amplitude: ');
fs = input('Enter sampling frequency: ');
ts = input('Enter time length [s]: ');
tiv =1/fs;
t = 0:tiv:ts-tiv;
carrier_signal = Ac*sin(2*pi*fc*t);
y = m.*carrier_signal; % according to definition

subplot(3,1,1),plot(t, m, 'r'), grid(), title('Message signal');
subplot(3,1,2), plot(t, carrier_signal, 'b'), grid(), title('Carrier Signal');
subplot(3,1,3), plot(t,y, 'g'), grid(), title('DSBSC');

end