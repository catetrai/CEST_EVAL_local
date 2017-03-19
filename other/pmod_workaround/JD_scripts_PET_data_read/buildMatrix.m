%BUILDMATRIX opens a Siemens Inveon (R) image file 
%
%   USAGE: 
%       buildMatrix(imageFile, headerInformation, ...
%                   downsample, subSection, saveFileName);
%           Will process the image given in 'imageFile' using the header-
%           information from headerReader.m in 'headerInformation'. 
%           Optional parameters:
%           - Downsample:   will load the data with a downsample factor. 
%                           i.e., 2, 4, 8, 16 (default: 1)
%                  CAUTION: Downsampling will skip voxels, it will NOT
%                           interpolate! This, to reduce loading time
%                           and reduce memory usage. For accurate
%                           downsampling: use the full matrix, and
%                           downsample manually afterwards.
%           - Subsection:   cut out one piece of the matrix 
%                           [startX,Y,Z,T; endX,Y,Z,T] (default: [1,1,1,1;
%                           Inf,Inf,Inf,Inf], i.e., everything)
%           - saveFileName: Do not load the data into memory but write as a
%                           file, as 'saveFileName'. Output will be the
%                           size of the matrix.
%
%           Output:
%           - imageMatrix:  the voxel matrix
%           - CTtestImage:  the dark and light scan, in case of a .cat file
%
% version 2016.10.13
% Last update: Bug fix, in the subsections.
%
% J.A. Disselhorst, 2009-2016
% University of Twente, Enschede (NL)
% Radboud University Medical Center, Nijmegen (NL)
% Werner Siemens Imaging Center, Tuebingen (DE)
%
% Disclaimer:
% THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT WARRANTY OF ANY
% KIND, EITHER EXPRESSED OR IMPLIED AND IS TO BE USED AT YOUR OWN RISK 

function [imageMatrix, testImages] = buildMatrix(imageFile,headerInformation,varargin)

    %% Input parsing -----------------------
    if nargin<2
        error('invalid number of arguments');
    else
        scanInfo = parseInputs(varargin);
    end
    N = size(scanInfo.subSection,2);
    if N<4
        scanInfo.subSection(1,N+1:4) = 1; 
        scanInfo.subSection(2,N+1:4) = Inf; 
    end
    testImages = 0; % default out.
    
    %% Determine the filetype and obtain file information --------------
    try 
        headerInformation.General.modality;
    catch 
        headerInformation.General.modality = -1;
    end
    scanInfo = obtainFileInformation(headerInformation,imageFile,scanInfo);
    if nargin<=2 % No downSample factor is given by the user, check required memory
        scanInfo = checkMemory(scanInfo);
    end
    
    %% Load the data
    switch scanInfo.typeOfData
        case 'Image'
            imageMatrix = loadImageFile(imageFile,scanInfo);
        case '3D Sinogram'
            if nargin>2
                warning('Downsampling and subsections are not available for 3D sinograms!');
            end
            imageMatrix = load3DSinogram(imageFile,scanInfo);
        case 'CatFile'
            if nargin>2
                warning('Downsampling and subsections are not available for 3D sinograms!');
            end
            [imageMatrix,testImages] = loadCatFile(imageFile,scanInfo);
    end
    
    %% Functions
    function scanInfo = obtainFileInformation(headerInformation,imageFile,scanInfo)
        data_type = headerInformation.General.data_type;
        if data_type>=5, scanInfo.endian = 'b'; else scanInfo.endian = 'l'; end;
        DataTypesFRead = {'int8=>int8','int16=>int16','int32=>int32','float32','float32','int16=>int16','int32=>int32'};
        DataTypesVar = {'int8',        'int16',       'int32',       'single', 'single', 'int16',       'int32'};
        DataSizes = [1 2 4 4 4 2 4];
        scanInfo.bytesize = DataSizes(data_type);
        scanInfo.datatypefread = DataTypesFRead{data_type};
        scanInfo.datatypevar = DataTypesVar{data_type};
        scanInfo.x = headerInformation.General.x_dimension;
        scanInfo.y = headerInformation.General.y_dimension;
        scanInfo.z = headerInformation.General.z_dimension;
        scanInfo.t = headerInformation.General.total_frames;
        scanInfo.scaleFactors = ones(1,scanInfo.t);
        for ii = 1:scanInfo.t
            try scaleFactor = eval(['headerInformation.frame_' num2str(ii-1) '.scale_factor;']); 
                if ischar(scaleFactor), scaleFactor = str2double(scaleFactor); end
                scanInfo.scaleFactors(ii) = scaleFactor;
            end
        end
        try scanInfo.calibrationFactor = headerInformation.General.calibration_factor;
        catch, scanInfo.calibrationFactor = 1; end;
        if ischar(scanInfo.calibrationFactor), scanInfo.calibrationFactor = str2double(scanInfo.calibrationFactor); end
        try scanInfo.isotopeBranchingFraction = headerInformation.General.isotope_branching_fraction;
        catch, scanInfo.isotopeBranchingFraction = 1; end
        try scanInfo.w = headerInformation.General.w_dimension;
        catch, scanInfo.w = 1; end
        
        fID  = fopen(imageFile,'r',scanInfo.endian);
        fseek(fID, 0, 'eof');
        scanInfo.filesize = ftell(fID);
        fseek(fID, 0, 'bof');
        fclose(fID);
        scanInfo.typeOfData = 'Image';
        switch headerInformation.General.modality
            case -1 % Unknown
                error('Unknown scanner modality. Aborting');
            case 0  % PET
                if scanInfo.w > 1  % 3D sinogram
                    scanInfo.typeOfData = '3D sinogram';
                end
            case 1  % CT
                if headerInformation.General.file_type == 14; % raw projection data
                    try 
                        scanInfo.z = headerInformation.General.number_of_projections;
                    catch
                        scanInfo.z = (scanInfo.filesize-headerInformation.General.ct_header_size)/scanInfo.bytesize/scanInfo.x/scanInfo.y;
                        fprintf('(Is it fluoroscopy imaging? - unkwown end results (total time frames %1.0f ??))\n',scanInfo.z);
                    end
                    scanInfo.ct_header_size = headerInformation.General.ct_header_size;
                    scanInfo.typeOfData = 'CatFile';
                end
        end
        
    end

    function scanInfo = checkMemory(scanInfo)
        if strcmpi(scanInfo.typeOfData,'Image')
            MatrixSize = scanInfo.x*scanInfo.y*scanInfo.z*scanInfo.t*scanInfo.bytesize;
            if MatrixSize>4E9
                fprintf('Image matrix will be large (%1.0f voxels)\n',MatrixSize)
                reply = questdlg(sprintf('Large image matrix (%1.0f voxels)\n',MatrixSize),'Memory','Ignore','Abort','Downsample', 'Abort');
                if isempty(reply) || strcmpi(reply,'abort')
                    ExitInveon;
                elseif strcmpi(reply,'downsample')
                    while MatrixSize>4E8
                        downFactor = downFactor*2;
                        MatrixSize = scanInfo.x*scanInfo.y*scanInfo.z*scanInfo.t*scanInfo.bytesize;
                    end
                end
            end
        else
            if scanInfo.downFactor>1
                warning('Down sampling only supported in images');
            end
        end
    end

    function imageMatrix = loadImageFile(imageFile,scanInfo)
        dims = [1,1,1,1; scanInfo.x, scanInfo.y, scanInfo.z, scanInfo.t];
        scanInfo.subSection = [max([dims(1,:); scanInfo.subSection(1,:)])-1; min([dims(2,:); scanInfo.subSection(2,:)])];
        scanInfo.subSection(3,:) = diff(scanInfo.subSection);
        if any(scanInfo.subSection(3,:)<1)
            error('Incorrect subsection');
        end
        
        ZSIZE = length(scanInfo.downFactor+scanInfo.subSection(1,3):scanInfo.downFactor:scanInfo.subSection(2,3));  
        matrixSize = [length(scanInfo.downFactor+scanInfo.subSection(1,1):scanInfo.downFactor:scanInfo.subSection(2,1)), ...
                      length(scanInfo.downFactor+scanInfo.subSection(1,2):scanInfo.downFactor:scanInfo.subSection(2,2)), ...
                      ZSIZE, scanInfo.subSection(3,4)];
                  
        % Display Information:
        fprintf('Loading file ''%s''\n',imageFile);
        fprintf('Size: %1.0f x %1.0f x %1.0f, in %1.0f frames\n',scanInfo.x,scanInfo.y,scanInfo.z, scanInfo.t);
        if scanInfo.downFactor>1
            fprintf('Downsampling factor: %1.0f\n',scanInfo.downFactor)
        end
        fprintf('Final matrix size: %1.0f x %1.0f x %1.0f x %1.0f\n',matrixSize);

        fID  = fopen(imageFile,'r',scanInfo.endian);
        if scanInfo.saveAs
            savefID = fopen(scanInfo.saveAs,'w+',scanInfo.endian);
        else
            imageMatrix = zeros(matrixSize, scanInfo.datatypevar);
        end
        
        for ii = scanInfo.subSection(1,4)+1:scanInfo.subSection(2,4)
            fprintf('Loading frame %1.0f: ....',ii);
            framebegin = -(scanInfo.x*scanInfo.y*scanInfo.z*(scanInfo.t-ii+1))*scanInfo.bytesize;
            slicebegin = scanInfo.x*scanInfo.y*scanInfo.subSection(1,3)*scanInfo.bytesize;
            fseek(fID, framebegin+slicebegin, 'eof'); %%%%%%%%%%%
            
            for j = 1:ZSIZE
                fseek(fID,scanInfo.x*scanInfo.y*(scanInfo.downFactor-1)*scanInfo.bytesize,'cof'); % Skip the downsample part
                tempImages = fread(fID,scanInfo.x*scanInfo.y,scanInfo.datatypefread);
                tempImages = reshape(tempImages,[scanInfo.x,scanInfo.y]);
                tempImages = tempImages(scanInfo.downFactor+scanInfo.subSection(1,1):scanInfo.downFactor:scanInfo.subSection(2,1),scanInfo.downFactor+scanInfo.subSection(1,2):scanInfo.downFactor:scanInfo.subSection(2,2));
                if scanInfo.saveAs
                    fwrite(savefID,tempImages,scanInfo.datatypevar);
                else
                    imageMatrix(:,:,j,ii-scanInfo.subSection(1,4)) = tempImages .* (scanInfo.calibrationFactor * scanInfo.scaleFactors(ii) / scanInfo.isotopeBranchingFraction);
                end
                fprintf('\b\b\b\b%3.0f%%',j/ZSIZE*100);
            end
            fprintf('\b\b\b\bComplete\n');
        end
        fclose(fID);
        if scanInfo.saveAs
            fclose(savefID);
            imageMatrix = matrixSize;
        end
    end

    function [imageMatrix,testImages] = loadCatFile(imageFile,scanInfo)
        fprintf('Loading CT .cat file....');
        fID  = fopen(imageFile,'r',scanInfo.endian);
        fseek(fID, scanInfo.ct_header_size, 'bof');
        testImages = fread(fID,scanInfo.x*scanInfo.y*2,scanInfo.datatypefread);
        testImages = reshape(testImages,[scanInfo.x,scanInfo.y,2]);
        
        imageMatrix = zeros(scanInfo.x, scanInfo.y, scanInfo.z, scanInfo.t, scanInfo.datatypevar);
        for ii = 1:scanInfo.t
            fprintf('\nLoading frame %1.0f: ....',ii);
            fseek(fID, -(scanInfo.x*scanInfo.y*scanInfo.z)*(scanInfo.t-ii+1)*scanInfo.bytesize, 'eof');
            for jj = 1:scanInfo.z
                tempImages = fread(fID,scanInfo.x*scanInfo.y,scanInfo.datatypefread);
                if ~isempty(tempImages)
                    tempImages = reshape(tempImages,[scanInfo.x,scanInfo.y]);
                    imageMatrix(:,:,jj,ii) = tempImages;
                end
                fprintf('\b\b\b\b%3.0f%%',jj/scanInfo.z*100);
            end
        end
        fprintf('\b\b\b\bComplete\n');
        fclose(fID);
    end

    function imageMatrix = load3DSinogram(imageFile,headerInformation,scanInfo)
        fprintf('Loading 3D sinogram....\n');
        fID  = fopen(imageFile,'r',scanInfo.endian);
        
        delta_elements = headerInformation.General.delta_elements;
        delta_sum = cumsum(delta_elements); 
        imageMatrix = zeros(scanInfo.x, scanInfo.y, delta_sum(end), scanInfo.t);
        
        for ii = 1:scanInfo.t
            fprintf('Loading frame %1.0f: ....',ii);
            fseek(fID, -(scanInfo.x*scanInfo.y*delta_sum(end)*(scanInfo.t-ii+1))*bytesize, 'eof');
            tempImages = fread(fID,scanInfo.x*scanInfo.y*delta_sum(end),scanInfo.datatypefread);
            imageMatrix(:,:,:,ii-subSection(1,4)) = reshape(tempImages,[scanInfo.x,scanInfo.y,delta_sum(end)]) .* (calibrationFactor * scaleFactor / isotopeBranchingFraction);
        end
        imageMatrix = Process3DSinogram(imageMatrix,delta_elements,scanInfo.w,currentFrame);
        
        fclose(fID);
    end

    function sino = Process3DSinogram(imageMatrix,delta_elements,wdim)
        fprintf('Processing 3D sinogram ...');
        [~,~,~,t] = size(imageMatrix);
        delta_sum = cumsum(delta_elements);
        sino = struct('data',[]);
        for ii = 1:t
            sino = setfield(sino, {ii,1}, 'data', squeeze(imageMatrix(:,:,1:delta_sum(1,2),ii)));
            sino = setfield(sino, {ii,1}, 'w', 0);
            for w = 2:wdim     %w_dimension;
                tempData = imageMatrix(:,:,delta_sum(w-1,2)+1:delta_sum(w,2),ii);
                sino = setfield(sino, {ii,w*2-1}, 'data', squeeze(tempData(:,:,end/2+1:end)));
                sino = setfield(sino, {ii,w*2-2}, 'data', squeeze(tempData(:,:,1:end/2)));
                sino = setfield(sino, {ii,w*2-1}, 'w', w-1);
                sino = setfield(sino, {ii,w*2-2}, 'w', -w+1);
            end
        end
        fprintf('\b\b\b\b: Complete\n');
    end

    function params = parseInputs(input)
        p = inputParser;
        p.addOptional('downFactor',1,@(x) isscalar(x) & round(x)==x & x>0);
        p.addOptional('subSection',[1 1 1 1; Inf Inf Inf Inf], @(x) all(size(x)==[2,4]) | all(size(x)==[2,3]));
        p.addOptional('saveAs',0,@ischar)
        p.parse(input{:});
        params = p.Results;
    end

    function ExitInveon()
        error('InveonFile:AbortedByUser','Aborted by user');
    end
end