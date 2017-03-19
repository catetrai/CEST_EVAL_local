function y = predict_pH(mdl,x)
% ** function y = predict_pH(mdl,x)
%
% 'mdl' is a fit object specifying the calibration linear model (pH ~ CESTratio).
% 'x' is the CESTratio map.
% 'y' is the predicted pH map.

imsize = size(x);
x = x(:);
b = mdl.Coefficients.Estimate;
y = [ones(numel(x),1) x] * b;
y(y<5 | y>8) = NaN;
y = reshape(y,imsize);