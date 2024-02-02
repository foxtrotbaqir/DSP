% An Elementary Single-Carrier M-QAM-based Digital Communication Chain
%% Basic system parameters
R = 10e6;
T = 1/R;
r = 4;
Fs = r/T;
Ts = 1/Fs;

%% Basic Transmitter Processing
N_symbols = 10000;
M = 4;

qam_axis = -sqrt(M)+1:2:sqrt(M)-1;
alphabet = qammod(0:sqrt(M)-1, M, 'unitaveragepower', true);
symbols = alphabet(randi(length(alphabet), 1, N_symbols));

% Plot the symbols
figure;
plot(symbols, 'ro', 'MarkerFaceColor', 'r');
axis equal;
xlabel('Real part');
ylabel('Imaginary part');
title('Transmitted symbols');
grid on;

%% Transmitter pulse shaping filtering
N_symbols_per_pulse = 40;
alpha = 0.20;

gt = rcosdesign(alpha, N_symbols_per_pulse, r, 'sqrt');
figure;
plot(-N_symbols_per_pulse*r/2*Ts:Ts:N_symbols_per_pulse*r/2*Ts, gt, 'b');
xlabel('Time [s]');
ylabel('Value');
title('Transmit RRC filter (pulse shape)');
grid on;

symbols_upsampled = upsample(symbols, r);
xt = filter(gt, 1, symbols_upsampled);

NFFT = 16384;
f = -Fs/2:1/(NFFT*Ts):Fs/2-1/(NFFT*Ts);

% Calculating and plotting the spectrum, in absolute scale and in dB
figure;
subplot(2,1,1);
plot(f/1e6, fftshift(abs(fft(xt, NFFT))));
xlabel('Frequency [MHz]');
ylabel('Relative amplitude');
title('TX signal amplitude spectrum');
grid on;

subplot(2,1,2);
plot(f/1e6, fftshift(20*log10(abs(fft(xt, NFFT)))));
xlabel('Frequency [MHz]');
ylabel('Relative amplitude [dB]');
title('TX signal amplitude spectrum');
grid on;
