clc;clear;
%{
    Open for educational purposes
    COMM.SYS.400 Digital Communication Matlab Project
    Baqir Kazmi
    %}

% TASK 2.1 TRANSMITTER
symbol_rate = 40e6;             % Rs = 40Mhz;
M = 8;                          % Modulation order = 8
roll_off = 0.1;                 % Roll-off factor alpha
oversample = 2;                 % Oversampling factor
Fs = oversample*symbol_rate;    % Sampling Freq = (Oversampling_factor)*(Symbol Rate)
Ts = 1/Fs;                      % Sampling Time interval
Nd = 20;                        % Duration of RRC Filter in symbols

org_bitstream = randi([0, 1], 1, (log2(M)*10^4)); % Random Bits generation
bitstream = org_bitstream;
%% The follwoing section implements Error control. To run without error coding comment out this section and Error decoder section from line 233 to 258.

%ERROR CONTROL CODING, HAMMING(7,4)
k = 4;                                      % number of source bits
n = 7;
ind = 1;                                    % variable to track for loop iterations
Par_arry = [1 0 1; 1 1 1; 1 1 0; 0 1 1];    % Parity matrix
G = [eye(k), Par_arry];                     % Generator matrix 
for i= 1:k:length(bitstream)
    cwdword(:,ind) = mod(bitstream(i:i+k-1)*G,2);   %Generating code words
    ind = ind+1;
end
bitstream = (cwdword(:))';

%%
bitstream = reshape(bitstream,log2(M),[]);          
%{ 
 The above line converts the singular bit stream into a log2(M)xY matix such that the each coloum represents a symbol and
the total number of coloums is equal to the total number of symbols 
%}
bit_dec = 2.^(log2(M)-1:-1:0);                      % Creates any array of powers of 2 to convert binary into decimal
symbol_stream = bit_dec*bitstream;                  % Converts the corresponding bits into their respective symbols indexis
Tx_signal = qammod(symbol_stream, M, pi/M,'gray');  % Applies PSK-modulation of order 'M' with phaseoffset = pi/8 and Gray mapping 
scatterplot(Tx_signal);grid on;                     % Plots the generated Tx symbols
title('16-QAM Constillation Diagram');

up_Txsignal = zeros(1, oversample*size(symbol_stream,2));                   % Creates an array of zeros
up_Txsignal(1:oversample:oversample*size(symbol_stream,2)) = Tx_signal;     % Inserts zeros between each symbol to oversample the signal
RRC = rcosdesign(roll_off, Nd, oversample, 'sqrt');                         % Generates an RRC filter with roll-off factor, Nd, oversample.
Filt_Tx = filter(RRC,1,up_Txsignal);                                        % Applying RRC-filter on usampled Tx signal
Filt_Tx  = Filt_Tx(1+(length(RRC)-1)/2:end);                                % Correcting filter delay

% Below lines of code is used to plot the Transmiter RRC-Filter response in time domain 
figure
freq_axis = Nd*oversample/2*Ts;
plot((-freq_axis:Ts:freq_axis).*(10^6),RRC,'k');hold on;
stem((-freq_axis:Ts:freq_axis).*(10^6),RRC,'ro');grid on;
xlabel('time [microsec]');ylabel('Amplitude');
title('Transmiter RRC filter Response');
legend('Pulse shape','Ideal symbol-sampling locations');

figure()
% Below function plots the Amplitude Response using FFT size of 2^14.
% (Function defined at the end) 
ampl_res(Filt_Tx, Fs, Ts);
title('Filtered TX signal');grid on;ylim([0 500]);

%%

%TASK 2.2
% Using Model 4 as our complex-valued channel impulse response.
h = [-0.4781 - 0.4781i, -0.3339 + 0.8264i, 0.0278 - 0.5302i, 0.2233 - 0.3868i, 0.0221 + 0.3155i];     % Model 4
L1 = 0;                             % Maximum tap is the first one, hence, L1 is 0 and L2 is 4. 
L2 = 4;
SNR = 0:2:25;

Propg_Tx = filter(h,1,Filt_Tx);     % Applying complex-valued channel to the Filtered transmitted signal.
Propg_Tx = Propg_Tx(L1+1:end);      % Removing the delay, if channel is anti-causal.
noise = (1/sqrt(2))*(randn(size(Propg_Tx))...
    + 1j*randn(size(Propg_Tx)));    % Complex white Gaussian random noise
P_sig = var(Propg_Tx);              % Signal power
P_noise = var(noise);               % Noise power
Prg_N_Tx = zeros(length(SNR), length(Propg_Tx));
noise_scaling_factor = sqrt(P_sig/P_noise./10.^(SNR./10)*(oversample/(1+roll_off))); % Scaling noise vector to get SNR 0 to 20 dB
for i = 1:length(noise_scaling_factor)
    Prg_N_Tx(i,:) = Propg_Tx +...
   (noise_scaling_factor(i)*noise);   %Transmitted signal + ISI + AWGN
end
for i = 1:length(SNR)
    figure(4)
    subplot(ceil(sqrt(length(SNR))), ceil(sqrt(length(SNR))), i)
    plot(real(Prg_N_Tx(i,:)), imag(Prg_N_Tx(i,:)),'.','LineWidth',3);grid on;  % Constillation plot of transmitted Tx symbols
    xlabel('In-Phase');ylabel('Quadrate Phase');
    title(sprintf('8-PSK(%d dB) + ISI + AWGN',SNR(i)));
end

%Plotting amplitude and phase response of the complex-valued channel:
figure(5);
[H,f]=freqz(h,1,-symbol_rate/2:symbol_rate/200:symbol_rate/2,symbol_rate); 
subplot(221)
plot(f/1e6,20*log10(abs(H)),'b'); grid on;
xlabel('Frequency [MHz]');ylabel('Amplitude response [dB]');
title('Frequency Responses');legend('ISI Channel');

subplot(222)
plot(f/1e6, phase(H)*(180/pi)); grid on;
xlabel('Frequency [MHz]');ylabel('Phase response [degree]');
title('Phase Responses');legend('ISI Channel');

figure
subplot(211)
ampl_res(Filt_Tx, Fs, Ts);legend('Amplitude response');
title('TX signal + RRC');grid on;ylim([0 500]);

subplot(212)
ampl_res(Propg_Tx, Fs, Ts); legend('Amplitude response');
title('TX signal + RRC+ ISI+ AWGN');grid on;ylim([0 500]);

%%
%TASK 2.3
Rx_signal = zeros(length(SNR), length(Prg_N_Tx));
rx_temp = 0;
for  i = 1:length(SNR)
    rx_temp = filter(RRC,1,Prg_N_Tx(i,:));                 % Applying RRC-filter on recieved Rx signal
    rx_temp = rx_temp(1+(length(RRC)-1)/2:end);            % Correcting filter delay
    rx_temp = downsample(rx_temp, oversample);             % Downsampling the recieved Rx signal
    Rx_signal(i,1:length(rx_temp)) = rx_temp;
end
Rx_signal(:, length(rx_temp)+1:end) = [];

for i = 1:length(SNR)
    figure(7)
    subplot(ceil(sqrt(length(SNR))), ceil(sqrt(length(SNR))), i)
    plot(real(Rx_signal(i,:)), imag(Rx_signal(i,:)),'.','LineWidth',3);grid on;  % Constillation plot of recieved Rx symbols
    xlabel('In-Phase');ylabel('Quadrate Phase');
    title(sprintf('Recieved 8-PSK(%d dB)',SNR(i)));
end


% Now implementing Least Mean Squares (LMS) algorithm-based equalizer for the channel equalization.
Equ_Rx_signal = zeros(length(SNR), length(Rx_signal));
N1 = 15;N2 = 15;
step_size = 0.001;                                   % step-size of the algorithm
LMS_coff = zeros(N1+N2+1,length(SNR));
for  j = 1:length(SNR)
    rx_temp = Rx_signal(j,:);
    c_LMS = zeros(N1+N2+1,1);                        % equalizer coefficients, initializations
    for i = N1+1:length(rx_temp)-N2
        rk = flipud(rx_temp(i-N1:i+N2).');           % Received signal vector
        Ek(i) = Tx_signal(i) - c_LMS.'*rk;           % Error signal, we assume a known symbol sequence
        c_LMS = c_LMS + step_size*Ek(i)*conj(rk);    % LMS update !
    end
    LMS_coff(:,j) = c_LMS;
    rx_temp = filter(c_LMS,1,Rx_signal(j,:));          % Applying channel estimation on recieved signal
    rx_temp = rx_temp(N1+1:end);                       % Correcting 
    Equ_Rx_signal(j,1:length(rx_temp)) = rx_temp;
end

figure(8);
hold on;grid on;
plot(abs(Ek)); title('Convergence behavior of the LMS-algorithm.');
ylabel('LMS error'); xlabel('Iteration index');

figure(9); hold on; stem(abs(conv(h,c_LMS)));grid on;
title('Effective impulse response (abs) of the equalized system ')

figure(10); 
[H,f]=freqz(h,1,-symbol_rate/2:symbol_rate/200:symbol_rate/2,symbol_rate); plot(f/1e6,20*log10(abs(H)),'b');
xlabel('Frequency [MHz]');ylabel('Amplitude response [dB]');

hold on; 
[H,f]=freqz(c_LMS,1,-symbol_rate/2:symbol_rate/200:symbol_rate/2,symbol_rate);
plot(f/1e6,20*log10(abs(H)),'r');

[H,f]=freqz(conv(c_LMS,h),1,-symbol_rate/2:symbol_rate/200:symbol_rate/2,symbol_rate);
plot(f/1e6,20*log10(abs(H)),'g');grid on;
legend('Channel','LMS Equalizer','Total Response (LMS)');

figure(5)
[H,f]=freqz(c_LMS,1,-symbol_rate/2:symbol_rate/200:symbol_rate/2,symbol_rate);
subplot(223)
plot(f/1e6,20*log10(abs(H)),'b'); grid on;
xlabel('Frequency [MHz]');ylabel('Amplitude response [dB]');
title('Frequency Responses');legend('ISI Channel');

subplot(224)
plot(f/1e6, phase(H)*(180/pi)); grid on;
xlabel('Frequency [MHz]');ylabel('Phase response [degree]');
title('Phase Responses');legend('LMS Equilizer');


Equ_Rx_signal(:, length(rx_temp)+1:end) = []; 
for i = 1:length(SNR)
    figure(11)
    subplot(ceil(sqrt(length(SNR))), ceil(sqrt(length(SNR))), i)
    plot(real(Equ_Rx_signal(i,:)), imag(Equ_Rx_signal(i,:)),'.','LineWidth',3);grid on;  % Constillation plot of Equilized Rx symbols
    xlabel('In-Phase');ylabel('Quadrate Phase');% legend(sprintf('Equilized 8-PSK- %d dB',SNR(i)));
    title(sprintf('Equilized 8-PSK- %d dB',SNR(i)));
end


% Plotting the original symbols, recieved symbols, and the equilized
% recieved symbols on the same plot for comparison.
for i = 1:length(SNR)
    figure(12);
    subplot(ceil(sqrt(length(SNR))), ceil(sqrt(length(SNR))), i)
    plot(real(Rx_signal(i,:)), imag(Rx_signal(i,:)),'y.','LineWidth',3);hold on;
    plot(real(Equ_Rx_signal(i,:)), imag(Equ_Rx_signal(i,:)), 'c.','LineWidth',3);hold on;
    plot(real(Tx_signal), imag(Tx_signal), 'kx','LineWidth',2);hold on;
    set(gca,'Color',[.7 .7 .7]);
    xlabel('In-Phase');ylabel('Quadrate Phase');title('Equilized 8-PSK Constillation Diagram');
end
legend('Recieved symbols','Equlized Recieved symbols','Original Transmitted symbols');grid on;
    
%%
demod_recieved = zeros(length(SNR), length(Equ_Rx_signal));
out_bitstream = zeros(length(SNR),length(demod_recieved)*3);              % Creating a variable to contain the output bitstream
for u = 1:length(SNR)
    rx_temp = [];
    rx_temp = pskdemod(Equ_Rx_signal(u,:),M,pi/M,'gray');       % Demodulating the recieved signal
    demod_recieved(u,:) = rx_temp;       % Demodulating the recieved signal
    for i = 1:length(demod_recieved)                              % Converting the symbol index into binary bitstream
        reminder = 0;
        for j = 1:log2(M)
            x = floor((demod_recieved(u,i)-reminder)/bit_dec(j));
            out_bitstream(u,((i-1)*3+j)) = x;
            reminder = reminder+bit_dec(j)*x;
        end
    end
    SNR_out(u,:) = out_bitstream(u,:); 
end
ber = zeros(1, length(SNR));

%% The follwoing section implements Error control. To run without error decoding comment out this section and Error coder section from line 20 to 33.
%ERROR CONTROL DECODING, HAMMING(7,4)

out_bitstream = []; 
for z = 1:length(SNR)
    bitstream2 = SNR_out(z,:);
    ind = 1;
    H = [Par_arry', eye(n-k)];
    error = mod(eye(n)*H',2);
    s = zeros(n-k, (floor(length(bitstream2)/n)*n-n)/n);
    recieved = zeros(k, (floor(length(bitstream2)/n)*n-n)/n);
    for i= 1:n:floor(length(bitstream2)/n)*n-n
        s(:,ind) = mod(bitstream2(i:i+n-1)*H',2);
        if sum(s(:,ind)) == 0
                recieved(:,ind) = bitstream2(i:i+k-1);
        else
            matching_rows = find(all(error == s(:,ind)', 2));
            bitstream2(i+matching_rows(1)-1) = ~bitstream2(i+matching_rows(1)-1);
            recieved(:,ind) = bitstream2(i:i+k-1);
        end
        ind = ind+1;
    end
    out_bitstream(z,:) = recieved(:)';
end


%%
for z = 1:length(SNR)
    err_vec = mod(org_bitstream(1:length(out_bitstream))+out_bitstream(z,:), 2); % Generating error vector
    err_vec = sum(abs(err_vec));                                    % Calculating total error
    ber(z) = err_vec/length(out_bitstream);                         % Calculating Bit-Error Rate 
    if ber(z) < 10^-8
        ber(z) = 10^-8; 
        % As 0 is undefined on a semi-log scale it cannot be plotted. Hence
        % BER is revalued as 10^-8 which will be the minimum BER value of
        % my semilog scale, just for illustration purposes.
        
    end
end

figure(13)
semilogy(SNR,ber, 'ko', 'LineWidth',2);grid on
title('Bit error rate');xlabel('SNR [dB]');ylabel('BER');
axis([min(SNR)-3 max(SNR)+3 10^-8 10]);

figure(14)
eyediagram(real(Filt_Tx), 8,  (1/symbol_rate));grid on;
title('Eye Diagram of Generated Signal');
xlabel('Time (s)');ylabel('Amplitude');


function ampl_res(signal, sampling_freq, sampling_time)     % Function to plot the Amplitude Response of a signal
    NFFT = 2^14;                                            % FFT size
    f = -sampling_freq/2:1/(NFFT*sampling_time):sampling_freq/2-1/(NFFT*sampling_time);  % frequency vector
    plot(f/1e6, fftshift(abs(fft(signal, NFFT))));
    xlabel('Frequency [MHz]');ylabel('Amplitude');
end
