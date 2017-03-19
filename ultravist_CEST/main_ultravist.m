% =========================== main_ultravist ==============================
% Script for pH mapping with Ultravist phantom (work in progress)
% CT 20170113

%% ------------------------- CEST-pH calibration --------------------------
%% Load phantom scans and fit z-spectra
clear directory; sequenceType = 'Mz'; READIMAGES = true;
load_CEST_bruker

Segment = make_Segment(M0_stack, 'ellipse', mean(M0_stack(M0_stack>0)).*[0.3]); 
dB0_stack_int = MINFIND_SPLINE_3D(Mz_stack, Segment, P);
Mz_CORR = B0_CORRECTION(Mz_stack, dB0_stack_int, P, Segment);
Z_corrExt = NORM_ZSTACK(Mz_CORR, M0_stack, P, Segment);
[Zlab, Zref, P, popt] = lorentzianfit_main(Z_corrExt, P, Segment, 'ultravist');

%% Calculate log10 ratio of MTR_Rex
delta = [4.2 5.6];
CESTratio = log10ratio(popt, P, Segment, delta);
figure; imagesc(CESTratio .* ~isnan(Segment)); colormap('jet'); colorbar;
title('log_{10}[MTR_{4.2ppm} / MTR_{5.6ppm}]')

%% Extract mean log10ratio from phantom ROIs
nPhan = 7;  % number of pH phantoms
dirfiles = what;
doesSegmentExist = ismember('Segment_ph.mat',dirfiles.mat); clear dirfiles;
if ~doesSegmentExist
    Segment_ph = nan(size(CESTratio,1), size(CESTratio,2), nPhan);
    for i=1:nPhan
        Segment_ph(:,:,i) = make_Segment(CESTratio, 'ellipse');
    end
else
    load('Segment_ph.mat');
end

CESTratio_ROIs = nan(1,nPhan);
CESTratio_ROIs_std = CESTratio_ROIs;
for i=1:nPhan
    roivals = CESTratio .* Segment_ph(:,:,i);
    CESTratio_ROIs(i) = nanmean(roivals(:));
    CESTratio_ROIs_std(i) = nanstd(roivals(:));
end
clear roivals

%% pH-CEST curve
ph =        [6.1   6.3   6.5   6.7   6.9   7.1   7.3];
% phantom#   7     6     5     4     3     2     1

CESTratio_ROIs = flip(CESTratio_ROIs);
figure;
plot(ph, CESTratio_ROIs, 'ok', 'markersize', 12,'linewidth',2)
% errorbar(ph, CESTratio_ROIs, CESTratio_ROIs_std, 'ko', 'markersize', 12,'linewidth',2)
set(gca,'xtick',ph);
xlabel('pH')
ylabel('log_{10}[MTR_{4.2ppm} / MTR_{5.6ppm}]')
title(sprintf('pH-CEST calibration at B1=%.1fuT', P.SEQ.B1))
set(gca,'fontsize',20)
hold on

% least-squares regression line
calibmdl = fitlm(ph,CESTratio_ROIs);
int_low = calibmdl.Coefficients.Estimate(1) - calibmdl.Coefficients.SE(1);
int_up = calibmdl.Coefficients.Estimate(1) + calibmdl.Coefficients.SE(1);
slope_low = calibmdl.Coefficients.Estimate(2) - calibmdl.Coefficients.SE(2);
slope_up = calibmdl.Coefficients.Estimate(2) + calibmdl.Coefficients.SE(2);

plot(ph, calibmdl.Fitted, '-r', 'linewidth', 2);
plot(ph, int_low + slope_low * ph, '--r', 'linewidth', 2);
plot(ph, int_up + slope_up * ph, '--r', 'linewidth', 2);

%% Linear regression fit
% calibration model (pH ~ CESTratio)
calibmdl = fitlm(CESTratio_ROIs, ph)
save('calibmdl.mat', 'calibmdl', '-v7.3');

%% --------------------------- In-vivo pH map -----------------------------
% Load and process in-vivo images using main 'analysis_CEST_bruker' script
% For multi-Lorentzian fitting, use:
[Zlab, Zref, P, popt] = lorentzianfit_main(Z_corrExt, P, Segment, 'invivo');

pos_4p2 = find_nearest(P.SEQ.w, 4.2);
pos_5p6 = find_nearest(P.SEQ.w, 5.6);
MTR_Rex = @(x,f) 1./Zlab(:,:,1,x) - 1./Zref.(f)(:,:,1,x);
CESTratiomap = log10(MTR_Rex(pos_4p2, 'ppm4p2') ./ MTR_Rex(pos_5p6, 'MT'));

pH_pred = predict_pH(calibmdl, CESTratiomap);
figure; imagesc(pH_pred); colorbar;
title('pH map')