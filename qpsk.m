function [qpsks] = qpsk(snr, T, f0, channelType)
    Fs = 2500e6; 
    N = 8192; 
    tler_28A = [1 1i 1 -1i 1 -1i -1 1i -1 1i 1 1i -1 1i 1 1i -1 1i 1 1i -1 1i -1 -1i 1 -1i 1 1i];
    p = T / 28.0;
    pn = fix(p * Fs); 
    cs = ones(28, pn);
    for j = 1:28
        if tler_28A(j) == -1
            cs(j, :) = cs(j, :) * -1;
        elseif tler_28A(j) == 1i
            cs(j, :) = cs(j, :) * 1i;
        elseif tler_28A(j) == -1i
            cs(j, :) = cs(j, :) * -1i;
        end
    end
    cs = [cs(1, :), cs(2, :), cs(3, :), cs(4, :), cs(5, :), cs(6, :), cs(7, :), ...
          cs(8, :), cs(9, :), cs(10, :), cs(11, :), cs(12, :), cs(13, :), cs(14, :), ...
          cs(15, :), cs(16, :), cs(17, :), cs(18, :), cs(19, :), cs(20, :), cs(21, :), ...
          cs(22, :), cs(23, :), cs(24, :), cs(25, :), cs(26, :), cs(27, :), cs(28, :)];
    t = (0:pn * 28 - 1) / Fs; 
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
    s = fft(s, N) / max(abs(fft(s, N)));
    qpsks.data = s;
    qpsks.label = 5;
end