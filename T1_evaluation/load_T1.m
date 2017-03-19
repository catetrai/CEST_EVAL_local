function T1_stack = load_T1(varargin)
% ** function T1_stack = load_T1(directory_TI_1, directory_TI_2, ... directory_TI_N)
%
% Inputs are the directory names of all Inversion Recovery images, one for
% each inversion time (TI).
% Output 'T1_stack' has dimensions (x,y,1,TI).
%
% CT 20170312

if nargin<1
    error('Must specify at least one image to load.')
end

nTI = numel(varargin);

% read metadata
size_xy = mineMetaDataFile([varargin{1} '/pdata/1/reco'], 'RECO_size');
size_zx = size_xy(1);
size_zy = size_xy(2);
% size_zz = mineMetaDataFile([directory_T1{1},'/method'], '##$PVM_NRepetitions')
size_zz = 1;

T1_stack = zeros(size_zx,size_zy,1,nTI);

id_read = waitbar(0,'Loading T1 mapping images...');
for i = 1:nTI
    waitbar(i/nTI);
    
    % read more metadata
    slope = mineMetaDataFile([varargin{i} '/pdata/1/reco'], 'RECO_map_slope');
    reco_word = mineMetaDataFile([varargin{i} '/pdata/1/reco'], 'RECO_wordtype');
    if strcmp(reco_word, '_32BIT_SGN_INT')
        reco_bit = 'uint32';
    elseif strcmp(reco_word, '_16BIT_SGN_INT')
        reco_bit = 'uint16';
    end
    
    fid = fopen([varargin{i} '/pdata/1/2dseq']);
    
    image_T1 = zeros(size_zx,size_zy,1,size_zz);
    for cont_z=1:size_zz
        for cont_i=1:size_zx
            for cont_j=1:size_zy
                image_T1(cont_i,cont_j,1,cont_z) = fread(fid,1,reco_bit);
            end
        end
    end
    
    image_T1 = image_T1/slope(1);   % scale image by correction factor    
    T1_stack(:,:,:,i) = mean(image_T1(:,:,1,:),4);
    
    fclose(fid);
end
close(id_read);