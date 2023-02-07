function p=powerspectrum( data, Fs )
    % Take the FFT of the data series and spit back meaningful power
    % spectrum formated data.
    T = 1/Fs;
    L = size( data , 1 );
    t = (0:L-1)*T;
    NFFT = 2^nextpow2(L);
    Y = fft(data,NFFT)/L;
    f = Fs/2*linspace(0,1,NFFT/2);
    power = 2*abs(Y(1:NFFT/2));
    p(1,:) = f;
    p(2,:) = power;
end