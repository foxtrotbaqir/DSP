function rx_demod = qpsk_demodulation(rx_dec)
% Demodulate the signal after decision
rx_demod = (rx_dec + 1)/2;
end