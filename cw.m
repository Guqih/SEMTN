function [cws] = cw(snr, T, f0, channelType)
    Fs = 2500e6; 
    N = 8192; 
    t = (0:1/Fs:T); 
    s = exp(1i * 2 * pi * f0 * t); 
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
    s = fft(s, N) / max(fft(s, N));
    cws.data = s;
    cws.label = 1;
end