function BER = norma(alg, Nt, Nr)
% MMSE algorithm
iter = 1000;
num = 120;
EbNo = 0:25;
idx = 1;

for SNR = EbNo
    errors = zeros(1, iter);
    for i=1:iter
        x = randn(1, num*Nt)>0;
        tx_bits = reshape(x, Nt, num);
        % Modulate the signal
        for a1 = 1:Nt
            tx_mod(a1,:) = qpsk_modulation(tx_bits(a1,:));
        end
        % Adding white noise as distortion
        H = (randn(Nr, Nt)+j*randn(Nr, Nt))/sqrt(2);
        [m,n] = size(tx_mod);
        spow = sum(sum(abs(tx_mod).^2))/(m*n);
        attn = 0.5 * spow/10.^(SNR/10);
        attn = sqrt(attn);
        noise = (randn(Nr, num)+j*randn(Nr, num))*attn;
        r = H * tx_mod + noise;
        % MMSE equalization
        if strcmpi(alg, 'MMSE')
            G = inv(H'*H+Nt/(10^(0.1*SNR))*eye(Nt))*H';
        end
        rx_equal = G*r;
        % De-modulation
        rx_dec = decision(rx_equal);
        for a2=1:Nt
            rx_demod(a2, :) = qpsk_demodulation(rx_dec(a2,:));
        end
        % Calculate the BER
        errors(i) = sum(sum(rx_demod~=tx_bits))/(Nt.*num);
    end
    BER(idx) = sum(errors)./iter;
    idx = idx+1;
end