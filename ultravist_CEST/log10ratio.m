function y = log10ratio(popt, P, Segment, delta)
% ** function y = log10ratio(popt, P, Segment, delta)
%
% Computes the log10 ratio of MTR_Rex at two specified frequency offsets
% (delta), for all pixels selected by the image segment.
%
% MTR_Rex is calculated as:
%       MTR_Rex(delta_i) = 1/Zlab(delta_i) - 1/Zref(delta_i)
%
% 'delta' is a vector of two frequency offset values for calculating the
%   ratio, with numerator followed by denominator.
%   e.g. for Ultravist, delta=[4.2 5.6]

[Zlab, Zref] = get_FIT_LABREF_modified(popt, P, Segment, 'ultravist', delta);
Zref_fields = regexprep({['ppm',num2str(delta(1))], ['ppm',num2str(delta(2))]},'\.','p');

MTR_Rex = @(x,f) 1./Zlab(:,:,1,x) - 1./Zref.(f)(:,:,1,x);
try
    y = log10(MTR_Rex(1, Zref_fields{1}) ./ MTR_Rex(2, Zref_fields{2}));
catch
    error('Zref field names do not correspond to the desired delta values.');
end