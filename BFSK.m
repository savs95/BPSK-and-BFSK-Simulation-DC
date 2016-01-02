% BER vs. SNR for FSK
%input = 1155;
enc7 = fec.bchenc(15,7);
enc11 = fec.bchenc(15,11);
decoder7 = fec.bchdec(enc7);
decoder11 = fec.bchdec(enc11);
fin_brr = [0,0,0,0,0,0,0,0];
fin_brr11 = [0,0,0,0,0,0,0,0];
fin_brr7 = [0,0,0,0,0,0,0,0];
%SNR=-20:3:10; %in dB
SNR=[-5 -2 0 2 5 8 12 15];
snr = 10.^(SNR/10); % actual ratio
N=1155; % number of bits
M=2; % number of levels in the signal.
Freq_sep=1600; % Hz
Nsamp = 8;
Fs =6400; % Hz
for i=1:5
    % Generate fskmod signal
    x = randi([0 M-1],N,1);
    x7 = encode(enc7,x);
    x11  = encode(enc11,x);
    y=fskmod(x,M,Freq_sep,Nsamp,Fs);
    y7=fskmod(x7,M,Freq_sep,Nsamp,Fs);
    y11=fskmod(x11,M,Freq_sep,Nsamp,Fs);
    signal_power=mean(abs(y).^2)*Nsamp;
    noise_power = (signal_power./snr);
    signal_power7=mean(abs(y7).^2)*Nsamp;
    noise_power7 = (signal_power7./snr);
    signal_power11=mean(abs(y11).^2)*Nsamp;
    noise_power11 = (signal_power11./snr);
    
    for k=1:length(SNR)
%         rcvda = awgn(y,SNR(k));
%         rcvd11a = awgn(y11,SNR(k));
%         rcvd7a = awgn(y7,SNR(k));
        
        noise = sqrt((noise_power(k)/2)).*(randn(1,N*Nsamp)+sqrt(-1)*randn(1,N*Nsamp));
        rcvd = y + noise';
        noise7 = sqrt((noise_power7(k)/2)).*(randn(1,N*Nsamp)+sqrt(-1)*randn(1,N*Nsamp));
        rcvd7 = y7 + noise7';
        noise11 = sqrt((noise_power11(k)/2)).*(randn(1,N*Nsamp)+sqrt(-1)*randn(1,N*Nsamp));
        rcvd11 = y11 + noise11';
        
        Z = fskdemod(rcvd,M,Freq_sep,Nsamp,Fs);
        Z11 = fskdemod(rcvd11,M,Freq_sep,Nsamp,Fs);
        Z7 = fskdemod(rcvd7,M,Freq_sep,Nsamp,Fs);
        Z7_d=decode(decoder7,Z7);
        Z11_d=decode(decoder11,Z11);
        errors(k) = sum(abs(x-Z));
        errors11(k) = sum(abs(x-Z11_d));
        errors7(k) = sum(abs(x-Z7_d));
    end
    fin_brr = fin_brr+(errors/N);
    fin_brr11 = fin_brr11+(errors11/N);
    fin_brr7 = fin_brr7+(errors7/N);
    
    
end
fin_brr = fin_brr/5;
fin_brr11 = fin_brr11/5;
fin_brr7 = fin_brr7/5;

semilogy(SNR,fin_brr,SNR,fin_brr11,SNR,fin_brr7);
legend('Without Error Correction','With BCH (15,11)','With BCH (15,7)');
title('BER vs. SNR plots for no coder and BCH (15,11), (15,7) coder over AWGN channel');
xlabel('SNR in dB')
ylabel('average BER in log scale')
grid