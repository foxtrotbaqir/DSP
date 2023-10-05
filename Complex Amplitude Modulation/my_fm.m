function y = my_fm(m,fc)
% where y = output
% x = message signal
% fc = carrier frequency
Ac = input('Enter carrier amplitude = ');
i = input('Enter Modulation Index = ');
fm = input('Enter message frequency = ');
fs = input('Enter sampling frequency: ');
ts = input('Enter time length [sec]: ');
tiv =1/fs;
t = 0:tiv:ts-tiv;
c = Ac*sin(2*pi*fc*t);
y = Ac*sin(2*pi*fc*t+i.*sin(2*pi*fm*t));
subplot(3,1,1), plot(t, m), xlabel('time'), ylabel('amplitude'), title('message');
subplot(3,1,2), plot(t, c), xlabel('time'), ylabel('amplitude'), title('carrier');
subplot(3,1,3), plot(t,y), xlabel('time'), ylabel('amplitude'), title('FM signal');
end