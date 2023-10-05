function y = my_hpf(x, fs, fcut)
% where 
% x = input
% n = order of filter
% fs = sampling frequency
% fcut = cutoff frequency 
% ftype = type of digital filter i.e. 'low', 'high', 'bandpass', 'stop', 'DC-0', 'DC-1' 
% y = output
n = 26;
ftype = 'high';
b = fir1(n,fcut,ftype);
y = filter(b,1,x);
% nfft = length(y);
% f = fs/2*linspace(0,1,nfft);
% plot(f,x)
% hold on
% plot(f,y)
% legend('Input Signal', 'Filtered Signal')

end

