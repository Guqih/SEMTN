function [qfsks] = qfsk(snr, T, f0, channelType)
    Fs = 2500e6; 
    N = 8192; 
    frank_16 = [0 0 0 0 0 1 2 3 0 2 4 6 0 3 6 9]; 
    p = T / 16.0;
    pn = fix(p * Fs); 
    cs = ones(16, pn);
    cst = (0:pn-1) / Fs; 
    df = (10 + (50 - 10) * rand(1)) * 1e6; 
    p0 = 0;
    for j = 1:16
        cs(j, :) = exp(1i * (p0 + 2 * pi * (f0 + df * mod(frank_16(j), 4)) * cst));
        p0 = p0 + 2 * pi * f0 * cst(pn); 
    end
    s = [cs(1,:), cs(2,:), cs(3,:), cs(4,:), cs(5,:), cs(6,:), cs(7,:), cs(8,:), ...
         cs(9,:), cs(10,:), cs(11,:), cs(12,:), cs(13,:), cs(14,:), cs(15,:), cs(16,:)];
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
    qfsks.data = s;
    qfsks.label = 7;
end