% note: parameters are selected after thorough calibration of functions and
% signals. they can be changed but every aspect where they are being used
% must be kept in mind.

% reading a sound file
samples = input('enter sampling rate of sound signal: ');
[s, s_rate] = audioread('sound_file.wav',[1, samples]);
figure;
plot(s)
NFFT = length(s);
S = fft(s,NFFT);
F = ((0:1/NFFT:1-1/NFFT)*s_rate).';
magnitudeS = abs(fftshift(S)); % Magnitude of the FFT
figure;
plot(F, magnitudeS), xlabel('Frequency axis'), ylabel('Magnitude'), title('Spectral Magnitude of Sound Signal');
%% Step 2 Filtering out unwanted frequencies
fcut = (1500/(s_rate/2));
sf = my_lpf(s,s_rate,fcut);
plot(sf);
SF = fft(sf,NFFT);
absSF = abs(fftshift(SF));
figure;
plot(F, absSF), xlabel('Frequency axis'), ylabel('Magnitude'), title('Spectral Magnitude of Filtered Signal');
% Note: In modulation and demodulation part, decomment or comment to use any mode
% of modulation or demodulation as required.

%% Step 3: modulating signal
% Modulation type: amplitude modulation
% sm1 = my_ampmod(sf, 500);
% sampling frquency must be higher atleast twice than fc and same as the no
% of samples taken in sound signal
% Modulation type: DSB-SC modulation
% sm2 = my_dsbsc(sf, 500);
% same parameter rules
% Modulation type: single side band
% sm3 = my_ssb(sf, 500);
% Modulation type: frequency modulation
 sm4 = my_fm(sf, 500);
% parameter rules same
%% Step 4: Demodulation
% demodulation type: coherent detection
% sd1 = my_coh(sm1, 500);
% demodulation type: envelope detection
%sd2 = my_env(sm4, 500);
% demodulation type: phase locked loop
sd3 = my_pll(sm4, 500);
% current setting is on mod type = fm and demod type = PLL

%% Step 5: Reconstructed Signal
sr = my_hpf(sd3,s_rate,fcut);
Sr = fftshift(fft(sr));
figure,
subplot(2,1,1), plot(sr), title('reconstructed signal');
subplot(2,1,2), plot(F, Sr), xlabel('Frequency [Hz]'), ylabel('Spatial Signal r'), title('Frequency domain [jw]');

%% --------------------------------------------------------------
% Section 6
% It is observed during performing various type of modulation and
% demodulation in the previous section that Phase locked Loop demodulation
% works the best in recovering lost signal with FM while amplitude
% modulation with envelope modulation has better result in recovering input
% signal than coherent detection.



