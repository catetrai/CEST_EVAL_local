function [imageMatrix,imageHeader,headerFile] = loadInveonFile(headerFile)
% LOADINVEONFILE opens an Inveon file
% 
% Usage: [imageMatrix,imageHeader,headerFile] = loadInveonFile(headerFile)
%
% Input: 
%         o headerFile: the location of a inveon header file [string]
%                       [optional, will ask for file when left blank]
%
% Output:
%         o imageMatrix: the image file [matrix of doubles]
%         o imageHeader: the header of the file [struct]
%         o headerFile: the location of a inveon header file [string]
%
% J.A. Disselhorst 2009-2014
% University of Twente, Enschede (NL)
% Radboud University Medical Center, Nijmegen (NL)
% Eberhard Karls University, Tuebingen (DE)
% Version 2014.12.01
%
% Disclaimer:
% THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT WARRANTY OF ANY
% KIND, EITHER EXPRESSED OR IMPLIED AND IS TO BE USED AT YOUR OWN RISK 

    if nargin<1 || isempty(headerFile)
        fprintf('Select header file: ...\n');
        PathName = getenv('JDisselhorstFolder');
        [FileName,PathName,FilterIndex] = uigetfile({'*.hdr','Header files (*.hdr)'},'Select Header file',PathName);
        headerFile = fullfile(PathName,FileName);
        if FilterIndex
            setenv('JDisselhorstFolder',PathName);
            clc
        else
            fprintf('Aborted\n');
            return
        end
    end
    imageHeader = headerReader(headerFile);
    if imageHeader.General.file_type > 1 % Not list-mode or unknown. 
        [imageMatrix, CTTestImages] = buildMatrix(headerFile(1:end-4),imageHeader);
        if ~CTTestImages, clear CTTestImages; 
        else
            k = menu('Use CT calibration files?', 'Yes', 'No');
            if k == 1
                imageMatrix = double(imageMatrix);
                CTTestImages = double(CTTestImages);
                subt = CTTestImages(:,:,1)./CTTestImages(:,:,2);
                for j = 1:size(imageMatrix,4)
                    for i = 1:size(imageMatrix,3)
                        imageMatrix(:,:,i,j) = imageMatrix(:,:,i,j)./CTTestImages(:,:,2) - subt;
                    end
                end
                clear subt k
            end
        end
    elseif imageHeader.General.file_type == 1 % List Mode
        fprintf('List-mode data currently not supported.\n');
        imageMatrix = [];
    end
end

