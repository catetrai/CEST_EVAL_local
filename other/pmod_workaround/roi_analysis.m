% initialize SPM functions
spm('defaults', 'fmri');
spm_jobman('initcfg');

% load DICOM of raw Mz image
mzhdr = dicominfo(fullfile(directory_Mz, '/pdata/1/dicom/MRIm01.dcm'));
hdr = mzhdr;

% change header
hdr.ImagesInAcquisition     = 1;
hdr.Filename                = '/Users/caterinatrainito/Documents/Caterina/Rat4_MTRrex_sl0.dcm';
hdr.PatientID               = 'Rat4_MTRrex_sl0';
hdr.StudyID                 = 'Rat4_MTRrex_sl0';
hdr.SmallestImagePixelValue = min(mtr_apt(:));
hdr.LargestImagePixelValue  = max(mtr_apt(:));
hdr.Slope                   = 6;

% change MTR image (get rid of NaNs)
mtr_apt(isnan(mtr_apt)) = 0;

% write new DICOM of CEST contrast image
dicomwrite(mtr_apt, '/Users/caterinatrainito/Documents/Caterina/Rat4_MTRrex_sl0.dcm', hdr);

% convert DICOM to Analyze using spm function
cd /Users/caterinatrainito/Documents/Caterina/MTRimages
dicom2img('/Users/caterinatrainito/Documents/Caterina/MTRimages/Rat4_MTRrex_sl0.dcm');

% read Analyze file with spm functions
petroi = spm_read_vols(spm_vol('/Users/caterinatrainito/Documents/Caterina/Rat4_PETroi.hdr'));


%% Reorient ROI for PET
clear matlabbatch

load('/Users/caterinatrainito/Documents/Caterina/transformmatrix/transfmat_attenmapbased_spm.mat');

% -------------------------------------------------------------------------
matlabbatch{1}.spm.util.reorient.srcfiles = ...
    {'/Users/caterinatrainito/Documents/Caterina/roi_analysis/Rat4/mask_healthy.img,1'};
matlabbatch{1}.spm.util.reorient.transform.transM = M;
matlabbatch{1}.spm.util.reorient.prefix = 'reorient_';
% -------------------------------------------------------------------------
spm_jobman('serial', matlabbatch);
clear matlabbatch

%% Interpolate ROI image to CEST and PET spaces
% Up-/down-sample ROI image with SPM batch ImCalc
% IMPORTANT: first image is the reference, second image is the ROI mask to
% be resampled
clear matlabbatch

roitype = {'tumor', 'healthy'};

ISPET=1;

for rt=1:2
% for rt=2
    if ISPET
        interptype = 1; % Trilinear (better for up-sampling to PET res)
    else
        interptype = 0; % NN (better for down-sampling to CEST res)
    end
    % -------------------------------------------------------------------------
    matlabbatch{1}.spm.util.imcalc.input = {
    '/Users/caterinatrainito/Documents/Caterina/dicoms_PET/sMaAm_For_14_06_Glio_dynFmiso_Rat4_10012017-0001-00178-000178.img,1'
    % '/Users/caterinatrainito/Documents/Caterina/MTRimages/sRat4_MTRrex_sl0-130001-00001-000001.img,1'
    % '/Users/caterinatrainito/Documents/Caterina/reorient_Rat4_singleSliceROI_ITK.img,1'
    % '/Users/caterinatrainito/Documents/Caterina/roi_analysis/Rat4/mask_healthy.img,1'
    sprintf('/Users/caterinatrainito/Documents/Caterina/roi_analysis/Rat4/reorient_mask_%s.img,1',roitype{rt})
                                            };
    matlabbatch{1}.spm.util.imcalc.output = sprintf('mask_%s_pet.img',roitype{rt});
    matlabbatch{1}.spm.util.imcalc.outdir = {'/Users/caterinatrainito/Documents/Caterina/roi_analysis/Rat4/'};
    matlabbatch{1}.spm.util.imcalc.expression = 'i2';
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = interptype;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    % -------------------------------------------------------------------------
    spm_jobman('serial', matlabbatch);
    clear matlabbatch
end

%% Sanity check CEST
cest = spm_read_vols(spm_vol('/Users/caterinatrainito/Documents/Caterina/MTRimages/sRat4_MTRrex_sl0-130001-00001-000001.hdr'));
cestroi = spm_read_vols(spm_vol('/Users/caterinatrainito/Documents/Caterina/roi_analysis/Rat4/mask_healthy_cest.hdr'));

if any(unique(cest))
    fprintf('\nCEST ROI has %d non-zero voxels. Values are between [%.2f, %.2f].\n\n', sum(cestroi(:)~=0), min(cestroi(:)), max(cestroi(:)));
else
    warning('CEST ROI has only null values!');
end

figure;
imagesc(cest .* ~cestroi);

%% Write CEST ROI values
cest(cest==0) = nan;          % make sure mask is not taking values outside of MTR brain-ROI
vals_healthy_cest = cest(logical(cestroi) & ~isnan(cest));
save('/Users/caterinatrainito/Documents/Caterina/roi_analysis/Rat4/vals_healthy_cest.mat', 'vals_healthy_cest');

%% Sanity check PET
pet = spm_read_vols(spm_vol('/Users/caterinatrainito/Documents/Caterina/dicoms_PET/sMaAm_For_14_06_Glio_dynFmiso_Rat4_10012017-0001-00178-000178.hdr'));
% petflip = spm_read_vols(spm_vol('/Users/caterinatrainito/Documents/Caterina/dicoms_PET/flippedPET178.hdr'));
petroi = spm_read_vols(spm_vol('/Users/caterinatrainito/Documents/Caterina/roi_analysis/Rat4/mask_healthy_pet.hdr'));

if any(unique(petroi))
    fprintf('\nPET ROI has %d non-zero voxels. Values are between [%.2f, %.2f].\n\n', sum(petroi(:)~=0), min(petroi(:)), max(petroi(:)));
else
    warning('PET ROI has only null values!');
end

nsl_pet = size(pet,3);
figure;
for sl = round(nsl_pet/3 : 2*nsl_pet/3)
    imagesc(pet(:,:,sl) .* ~flipud(petroi(:,:,sl)~=0));
    pause;
end

%% Write PET ROI values
petroi = petroi>0;          % PET mask must be binarized (because used trilinear interp)
petroi = flipud(petroi);    % L-R directions are flipped
vals_healthy_pet = pet(logical(petroi));
save('/Users/caterinatrainito/Documents/Caterina/roi_analysis/Rat4/vals_healthy_pet.mat', 'vals_healthy_pet');




%% Calculate slice number of T1 anatomical corresponding to CEST slice
mrdir=('/Users/caterinatrainito/Documents/Caterina/dicoms_MRanatom/sMaAm_FOR_14_06_Glio_dynFMISO_Rat4_10012017-530001-00001-000001.img');
mm2vox([0 0 0], mrdir)
mm2vox([0 0 1.5], mrdir)


%% LET'S DO THIS.

% load pet
[pet_img,pet_hdr,pet_path] = loadInveonFile('/Users/caterina/Desktop/WISC_MAC/PET_reconstr/MaAm_For_14_06_Glio_dynFmiso_Rat4_10012017/MaAm_For_14_06_Glio_dynFmiso_Rat4_10012017.img.hdr');
pet_img = squeeze(pet_img(:,:,:,2));

% apply transformation matrix derived from IRW or PMOD (manual
% coregistration)
tform = affine3d(tmat);
pet_trans = imwarp(pet_img,tform);
vis3d(pet_trans);

% load T1 post-gad
mr_img = nan(175,175,80);
for z=1:80
    mr_img(:,:,z) = dicomread(sprintf('/Users/caterina/Desktop/MaAm_FOR_14_06_Glio_dynFMISO_Rat4_10012017_1_1_20170110_100236/53/pdata/1/dicom/MRIm%02d.dcm',z));
end

% read pixel dump
vx_mr = textscan(fopen('/Users/caterina/Desktop/vx_mr.txt'), '%d %d %d');
vx_mr = double(cell2mat(vx_mr));
