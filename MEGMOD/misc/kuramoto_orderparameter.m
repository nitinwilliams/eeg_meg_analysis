function [R_mean,R_sigma]=kuramoto_orderparameter(D,mode)

num_samples=size(D,2);
Z=zeros(1,num_samples);

% Obtaining instantaneous phases
switch mode
    case 'hilbert'
    % D=unwrap(angle(hilbert(D')))';
    D=angle(hilbert(D'))';
    case 'morlet'
    % D=unwrap(angle(D'))';
    D=angle(D')';
end

for samp_idx=1:num_samples
    Z(samp_idx)=mean(exp(1i*D(:,samp_idx)));
end

R_mean=mean(abs(Z));
R_sigma=std(abs(Z));

end