function V = mock_hdr(P,h)
% ** function V = mock_hdr(P,h)
% Creates a mock header structure 'V' that can be used by SPM routines.
% Inputs:
%       'P' is the file name (full path) of the PET image that will be
%           written as '_NIFTI.img'. It is returned by JD's
%           'loadInveonFile' function as the variable 'headerFile'.
%       'h' is the header structure extracted using JD's 'loadInveonFile'
%           function.
%
% CT 20170127

% DATA TYPES (from the raw header text file)
% #   0 - Unknown data type
% #   1 - Byte (8-bits) data type
% #   3 - 4-byte integer - Intel style
% #   4 - 4-byte float - Intel style
% #   5 - 4-byte float - Sun style

%------------------------ MAPPING FROM spm_type ---------------------------
prec = {'uint8','int16','int32','float32','float64','int8','uint16','uint32'};
types= [    2      4       8       16        64       256    512      768];
%--------------------------------------------------------------------------

dtype = types(strcmp(prec,'float32'));
if h.General.data_type ~= 4
    warning('Data type is not single (assumed by default).');
end

V = struct('fname',   [P(1:end-8) '_NIFTI.img'],...
           'dim',     [h.General.x_dimension, h.General.y_dimension, h.General.z_dimension],...
           'dt',      [spm_type(dtype) 0],...
           'pinfo',   [1 0 0]',...
           'mat',     Vi(1).mat,...
           );
end