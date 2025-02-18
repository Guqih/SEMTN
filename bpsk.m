function [bpsks] = bpsk(snr, T, f0, channelType)
    Fs = 2500e6; 
    N = 8192; 
    barker_13 = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1]; 
    p = T / 13.0;
    pn = fix(p * Fs); 
    cs = ones(13, pn);
    for j = (1:13)
        if (barker_13(j) == -1)
            cs(j, :) = cs(j, :) * -1;
        end
    end
    cs = reshape(cs.', 1, []); 
    t = (0:pn * 13 - 1) / Fs; 
    s = exp(1i * 2 * pi * f0 * t) .* cs; 
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
    bpsks.data = s;
    bpsks.label = 4;
end