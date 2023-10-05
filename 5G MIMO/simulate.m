clear; clc;

BER1 = norma("MMSE", 2, 2);
BER2 = norma("MMSE", 2, 4);
BER3 = norma("MMSE", 4, 2);
BER4 = norma("MMSE", 4, 4);

EbNo = 0:25;

figure
semilogy(EbNo, BER1, 'y-x');
hold on
semilogy(EbNo, BER2, 'r-d');
hold on
semilogy(EbNo, BER3, 'g-*');
hold on
semilogy(EbNo, BER4, 'k-h');
hold on
xlabel('SnR [db]')
ylabel('BER')
grid on
title('Bit Error Rate with different size of mimo antennas using BPSK modulation')
legend('2x2 MIMO', '2x4 MIMO', '4x2 MIMO', '4x4 MIMO')