function y = my_ssb(m,fc)
% In radio communications, single-sideband modulation (SSB) or single-sideband 
% suppressed-carrier modulation (SSB-SC) is a type of modulation used to transmit
% information, such as an audio signal, by radio waves. A refinement of amplitude
% modulation, it uses transmitter power and bandwidth more efficiently. Amplitude
% modulation produces an output signal the bandwidth of which is twice the maximum
% frequency of the original baseband signal. Single-sideband modulation avoids this
% bandwidth increase, and the power wasted on a carrier, at the cost of increased
% device complexity and more difficult tuning at the receiver.
% where y = output signal
% m = message signal
% fc = carrier frequency
sb = input('Enter side of band, upper/lower: ');
Ac = input('Enter the value of carrier signal amplitude: ');
fs = input('Enter sampling frequency: ');
ts = input('Enter time length [sec]: ');

tiv =1/fs;
t = 0:tiv:ts-tiv;
mh = imag(hilbert(m));
if sb == 'upper'
    usb = m.*Ac.*cos(2*pi*fc*t) - mh.*Ac.*sin(2*pi*fc*t);
    y = usb;
    subplot(2,1,1), plot(t,y), xlabel('time'), ylabel('modulated signal');
    subplot(2,1,2), plot(abs(fftshift(fft(usb)))), xlabel('frequency'), ylabel('frequency response');
    
else
    lsb = m.*Ac.*cos(2*pi*fc*t) + mh.*Ac.*sin(2*pi*fc*t);
    y = lsb;
    subplot(2,1,1), plot(t,y), xlabel('time'), ylabel('modulated signal');
    subplot(2,1,2), plot(abs(fftshift(fft(lsb)))), xlabel('frequency'), ylabel('frequency response');
    
end

end