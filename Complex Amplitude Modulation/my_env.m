function z = my_env(y, fc)

% envelope detection
fc+2;
z = abs(hilbert(y));
plot(y, 'g'), title('modulated');
hold on;
plot (z , 'r'), title('envelope detected signal');
hold off;

end