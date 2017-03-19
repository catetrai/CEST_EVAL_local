function cest2dicom(cestImg, directory_Mz, dicomFilePath)
% ** function dicomhdr = cest2dicom(cestImg, directory_Mz, dicomFilePath)
%
% Writes a CEST contrast image created in MATLAB (e.g. MTR_Rex) into DICOM
% format. Also saves a .mat file of the scaling factor map used to convert
% double-precision values of the original CEST image to uint16 required by
% DICOM. Inputs:
%   'cestImg': the 2-D image array (double)
%   'directory_Mz': path of the raw Mz image (for DICOM header template)
%   'dicomFilePath': full path (including file name) of the DICOM file to
%       be written, e.g. '/Users/foo/Documents/exampleDicom.dcm'

% load DICOM header of Mz image
Mz_hdr = dicominfo(fullfile(directory_Mz, '/pdata/1/dicom/MRIm01.dcm'));
dicomhdr = Mz_hdr;
[~, PatientID] = fileparts(directory_Mz(1:max(strfind(directory_Mz,'/')-1)));

% change header
dicomhdr.ImagesInAcquisition     = 1;
dicomhdr.Filename                = dicomFilePath;
dicomhdr.PatientID               = PatientID;
dicomhdr.StudyID                 = PatientID;
dicomhdr.RescaleSlope            = 1;
dicomhdr = rmfield(dicomhdr, {'SmallestImagePixelValue', 'LargestImagePixelValue', 'FileModDate'});

% scale CEST image to uint16, and save the pixel-wise scaling factor map
cestImg_scaled = im2uint16(cestImg);
ScalingFactorMap = double(cestImg_scaled) ./ cestImg;
save([dicomFilePath(1:end-4) '_ScalingFactorMap.mat'], 'ScalingFactorMap', '-v7.3');

% write DICOM file
dicomwrite(cestImg_scaled, dicomhdr.Filename, dicomhdr);