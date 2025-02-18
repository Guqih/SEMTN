function [bfsks] = bfsk(snr, T, f0, channelType)
    Fs = 2500e6; 
    N = 8192; 
    barker_13 = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1]; 
    p = T / 13.0;
    pn = fix(p * Fs); 
    cs = ones(13, pn);
    cst = (0:pn-1) / Fs; 
    f1 = f0 + (10 + (50 - 10) * rand(1)) * 1e6; 
    p0 = 0; 
    for j = (1:13)
        if (barker_13(j) == 1)
            cs(j, :) = exp(1i * (p0 + 2 * pi * f0 * cst));
            p0 = p0 + 2 * pi * f0 * cst(pn);
        else
            cs(j, :) = exp(1i * (p0 + 2 * pi * f1 * cst));
            p0 = p0 + 2 * pi * f1 * cst(pn);
        end
    end
    s = reshape(cs.', 1, []); 
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
    bfsks.data = s;
    bfsks.label = 6;
end