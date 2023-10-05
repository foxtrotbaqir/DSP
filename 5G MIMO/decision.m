function rx_dec = decision(rx_equal)
% decides the type of signal
rx_denor = rx_equal;
rx_dec = sign(real(rx_denor));
end