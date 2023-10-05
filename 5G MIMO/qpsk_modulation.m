function tx_mod = qpsk_modulation(tx_bits)
% modulate the transmission signal on each antenna
table = [-1 1];
in = tx_bits;
tx_mod = table(in + 1);

end