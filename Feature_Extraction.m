clear;
clc;
disp('Generating');
N = 70000;
k = 0.5;
J = 8192;
ESR = 0;
%   channelType='Rayleigh + AWGN';
%   channelType='Rician + AWGN';
%   channelType='Nakagami + AWGN';
    channelType='AWGN';
for n=(0:0)
    T_min=0.5;T_max=3;
    f_min=200;f_max=220;
    snr_min=10;snr_max=10;
    for i=(1:(N/7))
        data1(i*7+1)=cw(snr_min+(snr_max-snr_min)*rand(1),(T_min+(T_max-T_min)*rand(1))*1e-6,(f_min+(f_max-f_min)*rand(1))*1e6,channelType);
        data1(i*7+2)=lfm(snr_min+(snr_max-snr_min)*rand(1),(T_min+(T_max-T_min)*rand(1))*1e-6,(f_min+(f_max-f_min)*rand(1))*1e6,channelType);
        data1(i*7+3)=nlfm(snr_min+(snr_max-snr_min)*rand(1),(T_min+(T_max-T_min)*rand(1))*1e-6,(f_min+(f_max-f_min)*rand(1))*1e6,channelType);
        data1(i*7+4)=bpsk(snr_min+(snr_max-snr_min)*rand(1),(T_min+(T_max-T_min)*rand(1))*1e-6,(f_min+(f_max-f_min)*rand(1))*1e6,channelType);
        data1(i*7+5)=qpsk(snr_min+(snr_max-snr_min)*rand(1),(T_min+(T_max-T_min)*rand(1))*1e-6,(f_min+(f_max-f_min)*rand(1))*1e6,channelType);
        data1(i*7+6)=bfsk(snr_min+(snr_max-snr_min)*rand(1),(T_min+(T_max-T_min)*rand(1))*1e-6,(f_min+(f_max-f_min)*rand(1))*1e6,channelType);
        data1(i*7+7)=qfsk(snr_min+(snr_max-snr_min)*rand(1),(T_min+(T_max-T_min)*rand(1))*1e-6,(f_min+(f_max-f_min)*rand(1))*1e6,channelType);
        x=['n=',num2str(n),'i=',num2str(i)];
        disp(x);
    end
end
feature = zeros(7, N); 
Labels = zeros(1, N);
for i = 1:N
    signal = data1(i+7).data;
    Mag = abs(signal); 
    feature(1, i) = std(Mag); 
    feature(2, i) = mean(Mag); 
    feature(3, i) = min(Mag); 
    Thr = k * mean(Mag); 
    ESR = sum(Mag > Thr) / J;
    feature(4, i) = ESR; 
    X = fft(signal); 
    phase_spectrum = angle(X); 
    phase_diff = diff(phase_spectrum); 
    feature(5,   i) = mean(phase_spectrum); 
    feature(6, i) = mean(phase_diff); 
    feature(7, i) = std(phase_diff); 
    Labels(i)=data1(i+7).label; 
end
disp('Finished');