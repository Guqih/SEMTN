function [nlfms] = nlfm(snr, T, f0, channelType)
Fs=2500e6;
N=8192;
B=(45+(55-45)*rand(1))*1e6;
syms f;
W(f)=0.54-0.46*cos(2*pi/B*f);
Tf=matlabFunction(vpa(double(vpa(T/int(W,-B/2,B/2)))*int(W,-B/2,f)));
f=(-B/2 : B/(T*Fs) : B/2);
[xData, yData] = prepareCurveData( Tf(f), f );
Ft= fit( xData, yData, 'smoothingspline');
t=(0 : 1/Fs : T);
Ft=Ft(t);
P=zeros(1,length(t));
for j=2:length(t)
    P(j)=P(j-1)+2*pi*Ft(j)*(1/Fs);
end
s=exp(1i*(2*pi*f0.*t+P));
    switch lower(channelType)
        case 'rayleigh + awgn'
            rayleighChan = comm.RayleighChannel('SampleRate', Fs, ...
                                                'PathDelays', [0, 5e-8], ...
                                                'AveragePathGains', [0, -8], ...
                                                'MaximumDopplerShift', 2);
            s = rayleighChan(s.'); 
            s = awgn(s, snr); 
        case 'rician + awgn'
            ricianChan = comm.RicianChannel('SampleRate', Fs, ...
                                            'PathDelays', [0, 1e-6], ...
                                            'AveragePathGains', [0, -8], ...
                                            'KFactor', 20, ...
                                            'MaximumDopplerShift', 4);
            s = ricianChan(s.'); 
            s = awgn(s, snr); 
        case 'nakagami + awgn'
            m = 2; 
            omega = 1; 
            h = sqrt(gamrnd(m, omega / m, size(s))); 
            s = s .* h;
            s = awgn(s, snr);
        case 'awgn'
            s = awgn(s, snr);
    end
    s = fft(s, N) / max(abs(fft(s, N)));
    nlfms.data = s;
    nlfms.label=3;
end