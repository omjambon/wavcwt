clear all
M = 4;          % Size of the signal constellation
k = log2(M);    % Number of bits per symbol

samp_freq = 40000000; %in Hz
fs=samp_freq;
signal_duration = 2;  %in seconds
t = 0:(1/(samp_freq)):signal_duration;

bs = randi([0 1],500,1); % Input signal
x = randi([0 M-1],50000,1); % Input signal
xxx = matfile('xxx.mat');
% bin_x = bi2de(reshape(x,k,length(x)/k).','left-msb');
bsg = bin2gray(x,'fsk',4);
modData = fskmod(bsg,M,4,4,samp_freq,'cont','gray');
frequency_sequence = randi([1 8],50000,1); % Input signal
frequency_sequence = [12;5;14];
% frequency_sequence = [4;7;1;3;5;7;5;2;7;4;8;3;5;6;4;8;2;3;7;1;5;2;3;7;6;3;7;8;5;2;1;4;6;2];



%% From binary sequence to gray mapping
bss = [0; 0; 1; 0; 1; 1; 1; 0; 0; 1];
% bin_x = bi2de(reshape(bss,k,length(bss)/k).','left-msb');
% bsg = bin2gray(bin_x,'fsk',4);

% degray_bsg = gray2bin(bsg,'fsk',4);
% debin_x = de2bi(degray_bsg,[],'left-msb');
% debin_x = reshape(debin_x',[],1);

%%
new_sig = [];
receiveData = [];

ber = [];
EbN0 = [15:35];
sig_dehop = [];
stock_sig_dehop = zeros(20,325);
receiveData = awgn(modData,-5,1);
% receiveData = awgn(modData,EbN0(k),1);
signal_duration = 2.5;  %in seconde %0.3
t = 0:(1/samp_freq):signal_duration;
for i = 1 : 2
    new_sig(i*64+1:i*64+64) = (modData(i*64+1:i*64+64)').*exp(1j*2*pi*frequency_sequence(i+1)*1000000*t(i*64+1:i*64+64));
    new_sig(3*64+1:3*64+64) = 0
    new_sig(4*64+1:4*64+64) = 0
% sig_dehop(i*100+1:i*100+100) = (new_sig(i*100+1:i*100+100)).*exp(-1j*2*pi*frequency_sequence(i+1)*1000000*t(i*100+1:i*100+100));
end

%% Ondelette
output1 = [real(new_sig);imag(new_sig)];
filename ='C:\Users\HP\Desktop\wavelet\test_python2.iq';
fid = fopen(filename,'w');
fwrite(fid,new_sig,'float');
% dlmwrite('C:\Users\charl\Desktop\textssss.txt',new_sig,'delimiter',' ')
% dlmwrite('C:\Users\charl\Desktop\textssss.txt',new_sig)

fclose(fid);
% spectrogram(new_sig,'yaxis')
cwt(new_sig,40000000)
[minf,maxf] = cwtfreqbounds(numel(new_sig),20e6);


%% AWGN WORK
snr = 25;
noisy_sig = awgn(new_sig,snr);
output_noisy_sig = [real(new_sig);imag(new_sig)];
filename ='C:\Users\HP\Desktop\wavelet\pasnoisy.iq';
fid_noisy = fopen(filename,'w');
fwrite(fid_noisy,output_noisy_sig,'double');