function [Z_stack, B1_input] = make_5D_B1_stack(protocol, dB0_stack_ext, Segment)
% ** function [Z_stack, B1_input] = make_5D_B1_stack(protocol, dB0_stack_ext, Segment)
%
% Creates 5D stack (x,y,z,w,B1) of z-spectra for each B1 amplitude.
%
% CT 20170111

sequenceType = 'Mz'; READIMAGES = false;
load_CEST_bruker
M0_stack = load_M0(protocol(ix.M0).Path);

nB1 = length(ix.Mz);
B1_input = nan(1,nB1);

for i = 1:nB1
    fprintf('Processing Z-spectrum for B1 = %.1f muT ...\n', protocol(ix.Mz(i)).B1);
    B1_input(i) = protocol(ix.Mz(i)).B1;
    directory_Mz = protocol(ix.Mz(i)).Path;
    Mz_stack = load_Mz(directory_Mz);
    P = wipread_modified(directory_Mz, protocol(ix.M0).Path);

    % calculate B0-corrected z-spectra using B0 map from WASABI
    Mz_CORR = B0_CORRECTION(Mz_stack, dB0_stack_ext, P, Segment);
    Z_corrExt = NORM_ZSTACK(Mz_CORR, M0_stack, P, Segment);

    Z_stack(:,:,:,:,i) = Z_corrExt;
end
fprintf('Done.\n')