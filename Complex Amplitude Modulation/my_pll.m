function z = my_pll(y,fc)


fv = input('enter oscillating frequency of VCO: ');
kv = input('enter gain of VCO: ');
fs = input('enter sampling frequency: ');
nf = input('enter number of samples in simulation: ');
n = input('enter number of order of filter: ');
% designing a lowpass filter
fcut = fc/(fs/2);
b = fir1(n,fcut);

ts = 1/fs;
t = [0:ts:(nf-1)*ts];

% VCO design
vco = zeros(1,nf);
phi = zeros(1,nf);
ref = y; % input signal

for i = 2:nf
    time = (n-2)*ts;
    error_mult(i) = ref(i)*vco(i-1);

    % passing into filter raw error
    for j=1:length(b)
        if i-j+1>=1
            error_array(j)=error_mult(i-j+1);
        else
            error_array(j)=0;
        end
    end
    error(i)=sum(error_array.*(b));
    phi(i) = phi(i-1)+2*pi*error(i)*kv*ts;
    vco(i) = sin(2*pi*fv*time+phi(i));
end
% plotting input and output signals
z = vco;
figure;
subplot(2,1,1), plot(t,ref,'b'), title('Plot of input signal'), xlabel('time [s]'), ylabel('amplitude');
subplot(2,1,2), plot(t,z,'r'), title('Plot of output signal'), xlabel('time [s]'), ylabel('amplitude');



end
