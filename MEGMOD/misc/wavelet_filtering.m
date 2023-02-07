function [D_wavfilt]=wavelet_filtering(D,m,fc,fs)

[M,N]=size(D);
D_wavfilt=zeros(M,N);

% Initialising values
LB=-1;
UB=1;
wavt=[LB:1/fs:UB];
s=m/(2*pi*fc);
A=1/sqrt(s*sqrt(pi));

% Generating Morlet wavelet
morwav=A*exp(2*pi*1i*fc*wavt).*exp((-wavt.^2)./(2*s^2));

half_wavlen = ceil(length(morwav)/2);
n_conv = length(morwav) + size(D,2) - 1;
fft_w = fft(morwav,n_conv);
fft_w_norm = fft_w./abs(max(fft_w));

% Multiplication in frequency domain

for idx=1:M
 fft_e = fft(D(idx,:),n_conv);
 D_wavfiltvec = ifft(fft_e.*fft_w_norm,n_conv);
 % D_wavfiltvec = ifft(fft_e.*fft_w,n_conv);
 D_wavfilt(idx,:) = D_wavfiltvec(half_wavlen:end-half_wavlen+1);
 % D_wavfilt(idx,:) = conv(D(idx,:),morwav,'same');
end

end