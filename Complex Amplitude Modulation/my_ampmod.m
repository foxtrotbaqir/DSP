function y = my_ampmod(m, fc)
% where m = message signal
% fc = carrier frequency must be really high
% y = output signal
% it is to be noted that time frame of both signals should be same
% otherwise an error will be generated.
fs = input('Enter sampling frequency: ');
Ac = input('Enter carrier amplitude: ');
ka = input('Enter amplitude sensitivity: ');
ts = input('Enter time length: ');
tiv = 1/fs; %time interval value
t = 0:tiv:ts-tiv;
c = Ac*cos(2*pi*fc*t);
y = (1+ka*m).*c; % modulation formula
figure;
subplot(3,1,1), plot(t,m, 'r'), xlabel('time'), ylabel('amplitude');
subplot(3,1,2), plot(t,c, 'g'), xlabel('time'), ylabel('amplitude');
subplot(3,1,3), plot(t, y, 'b'), xlabel('time'), ylabel('amplitude');
end

