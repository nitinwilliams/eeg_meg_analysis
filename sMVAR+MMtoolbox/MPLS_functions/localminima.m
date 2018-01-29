function indmin = localminima(x);
% indmin = localminima(x);

if (length(size(x)) > 2) || (numel(x) < 3)
    error('x is expected to be a column vector of three or more elements.');
end
npoints = length(x);
lpoints = zeros(npoints, 1);
fleq_itsleft = (x(2:end-1) < x(1:end-2));
fleq_itsright = (x(2:end-1) < x(3:end));
lpoints(2:end-1) = (fleq_itsleft & fleq_itsright);
indmin = find(lpoints);
return;