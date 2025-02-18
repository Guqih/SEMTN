function [lfms] = lfm(snr, T, f0, channelType)
    Fs = 2500e6; 
    N = 8192; 
    t = (0:1/Fs:T); 
    B = (45 + (55 - 45) * rand(1)) * 1e6; 
    K = rand(1) - 0.5; 
    if K > 0
        K = B / T; 
    else
        K = -B / T; 
    end
    s = exp(1i * (2 * pi * f0 * t + pi * K * t.^2)); 
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
    lfms.data = s;
    lfms.label = 2;
end