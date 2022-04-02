fid = fopen('C:\Users\HP\Desktop\wavelet\bpsk_in_noise.iq', 'r');  %for read
data = fread(fid);  %read two bytes at a time
Fs= 40e6

t = 0:1/Fs:1;
% tt = timetable(data(:),'RowTimes',seconds(t'));
% cwt(tt)

fb = cwtfilterbank(SignalLength=numel(data),SamplingFrequency=Fs,...
    FrequencyLimits=[10e6 20000000]);
% freqz(fb)
figure;
cwt(data,Fs)