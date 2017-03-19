%% CEST EVAL
% MODIFIED FROM CEST_EVAL toolbox:
%   Date: 2016/02/10 
%   Version for CEST-sources.de
%   Authors: Moritz Zaiss  - m.zaiss@dkfz.de , Johannes Windschuh 
%   Divison of Medical Physics in Radiology
%   German Cancer Research Center (DKFZ), Heidelberg Germany, http://www.dkfz.de/en/index.html
%   CEST sources - Copyright (C) 2016  Moritz Zaiss
%   **********************************
%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or(at your option) any later version.
%    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%    You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
%   **********************************


%% Add main CEST code folder and subfolders to path
code_directory = '~/Documents/MATLAB/CEST_EVAL-master-2';
cd(code_directory);
addpath(genpath('.'));

% prepare for loading data
clear directory

%% =================== FIELD INHOMOGENEITY CORRECTIONS ====================
%% ------------ (A) If WASABI sequence has *NOT* been acquired ------------
% (1) Load Mz and M0
sequenceType = 'Mz';
load_CEST_bruker

% (2) Define 2D Segment mask of ones and NaNs
Segment = make_Segment(M0_stack, 'free', mean(M0_stack(M0_stack>0)).*[0.3]);
%Segment = ones(size(M0_stack));

% (3) Calculate internal dB0 map
dB0_stack_int = MINFIND_SPLINE_3D(Mz_stack, Segment, P);

% (4) B0-correction of raw Mz image using internal B0-map
Mz_CORR = B0_CORRECTION(Mz_stack, dB0_stack_int, P, Segment);

% (5) Calculate B0-corrected Z-spectra for each pixel
Z_corrExt = NORM_ZSTACK(Mz_CORR, M0_stack, P, Segment);

%% IMGUI: look at raw Z-spectra
close(imgui); imgui

%% ------------- (B) If WASABI sequence *HAS* been acquired ---------------
% (1) Load WASABI stack (as Mz) and M0
sequenceType = 'WASABI';
load_CEST_bruker

% (2) Define 2D Segment mask of ones and NaNs
Segment = make_Segment(M0_stack, 'free', mean(M0_stack(M0_stack>0)).*[0.3]); 
% Segment = ones(size(M0_stack));

% (3) Obtain B0 map and B1 map using WASABI method
% Calculate WASABI spectra
Z_uncorr = NORM_ZSTACK(Mz_stack, M0_stack, P, Segment);

P.FIT.options   = [1E-04, 1E-15, 1E-10, 1E-4, 1E-06];
P.FIT.nIter     = 100;
P.FIT.modelnum  = 021011; % 021011 = WASAB1

[popt, P] = FIT_3D(Z_uncorr, P, Segment);

B1map = popt(:,:,1) / P.SEQ.B1;     % B1 map
dB0_stack_ext = popt(:,:,2);        % B0 map

%% See WASABI spectrum (choose map 'Z_uncorr', click on pixel) and B0/B1 maps:
close(imgui); imgui

%% (4) Create 5D stack (x,y,z,w,B1) of Z-spectra for different B1 values
%     Select all Mz scans acquired at different B1
[Z_stack, B1_input] = make_5D_B1_stack(protocol, dB0_stack_ext, Segment);

% (5) Specify the B1 value you want to reconstruct (i.e. the Mz acquisition
%     you will use for further analyses)
%     NOTE: will overwrite WASABI Mz image and parameter structure
sequenceType = 'Mz';
load_CEST_bruker
B1_output = protocol(ix.Mz).B1;
%B1_output = 1;  % in uT

% (6) Run B1 correction
Z_stack_corr = Z_B1_correction(Z_stack, B1map, B1_input, B1_output, Segment, 'linear');
Z_corrExt = Z_stack_corr(:,:,:,:,1);


%% ======================= MULTI-LORENTZIAN FITTING =======================

%For Ultravist phantom analysis, use function instead:
%[Zlab, Zref, P, popt] = lorentzianfit_main(Z_corrExt, P, Segment, 'ultravist')

P.FIT.options   = [1E-04, 1E-15, 1E-10, 1E-04, 1E-06];
P.FIT.nIter     = 50;
P.FIT.modelnum  = 5;     % number of Lorentzian pools (possible 1-5)
P.FIT.extopt    = 1;     % change parameters explicitly

% Lorentzian line Li defined by amplitude Ai, width Gi [ppm] and offset dwi[ppm]:
% 
%                   Li = Ai.*Gi^2/4 ./ (Gi^2/4+(dw-dwi).^2)
%
% 5-pool model: 1=water; 2=amide; 3=NOE; 4=MT; 5=amine
%
% const.Zi   A1   G1   dw1      A2     G2    dw2     A3    G3   dw3     A4    G4   dw4    A5    G5   dw5
lb = [ 0.5  0.02  0.3  -1       0.0    0.4   +3     0.0    1    -4.5    0.0   10   -4     0.0   1    1   ];
ub = [ 1    1     10   +1       0.2    3     +4     0.4    5    -2        1   100  -1     0.2   3.5  2.5 ];
p0 = [ 1    0.9   1.4   0       0.025  0.5   3.5    0.02   3    -3.5    0.1   25   -1     0.01  1.5  2.2 ];
P.FIT.lower_limit_fit = lb; P.FIT.upper_limit_fit = ub; P.FIT.start_fit = p0;

% Perform the fit pixelwise
[popt, P] = FIT_3D(Z_corrExt, P, Segment);
  
% Create Reference values Z_Ref=(Z_lab - Li)
[Zlab, Zref] = get_FIT_LABREF(popt, P, Segment);

%% IMGUI: look at fitted Z-spectra with individual Lorentzian components
close(imgui); imgui


%% ============================== T1 MAPPING ==============================
sequenceType = 'T1mapping';
load_CEST_bruker

% Information about fit
P_T1.SEQ.w         = [protocol(ix.T1mapping).TI]';      % inversion times
P_T1.FIT.options   = [1E-04, 1E-15, 1E-10, 1E-04, 1E-06];
P_T1.FIT.nIter     = 100;
P_T1.FIT.modelnum  = 031011;

% Starting parameters (optional)
%     T1         a      c
lb = [0         -5000   0       ];
ub = [5000      5000   5000    ];
p0 = [1000      -2000   1000    ];
P_T1.FIT.lower_limit_fit = lb; P_T1.FIT.upper_limit_fit = ub;
P_T1.FIT.start_fit = p0; StartValues=p0;

% Run T1 mapping
nROIs = 1;
[T1info, T1map, popt_T1] = T1eval_levmar(T1_stack,1,nROIs,P_T1,Segment,StartValues);


%% ====================== CEST CONTRAST CALCULATION =======================
% Change 'Amide' and relative chemical shift value (3.5) to whichever
% peak/pool you want to look at.

pos_Amide = find_nearest(P.SEQ.w, 3.5);     % find index of Amide peak

% MTR_LD: MTR asymmetry as linear difference
MTR_LD = @(x,f) Zref.(f)(:,:,1,x) - Zlab(:,:,1,x);
contrast_img = MTR_LD(pos_Amide, 'Amide');

% MTR_Rex: MTR asymmetry as difference of reciprocals
MTR_Rex = @(x,f) 1./Zlab(:,:,1,x) - 1./Zref.(f)(:,:,1,x);
contrast_img = MTR_Rex(pos_Amide, 'Amide');

% AREX: T1-corrected contrast
AREX = @(x,f) (1./Zlab(:,:,1,x) - 1./Zref.(f)(:,:,1,x)) ./ T1map;
contrast_img = AREX(pos_Amide, 'Amide');

figure; imagesc(contrast_img); colorbar


%% ====================== WRITE CEST CONTRAST IMAGE =======================
% Set full path (including file name .dcm) of DICOM file to be written
dicomFilePath = '';
cest2dicom(contrast_img, protocol(ix.Mz).Path, dicomFilePath);


