function dicom2img(files)
% function dicom2img(files)
%
% This function is to convert dicom files to Analyze image files using SPM
% script
%
% files (optional): list of dicom files
%
% Xu Cui
%

if(nargin < 1)
    files = spm_select('list', pwd, '\.dcm');
end

spm_get_defaults;
hdr = spm_dicom_headers(files);
spm_dicom_convert(hdr);
display('done!')
return;