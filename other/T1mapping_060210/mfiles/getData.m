% getData.m
%
% Generates the .mat from the DICOM images. 
%
% written by J. Barral, M. Etezadi-Amoli, E. Gudmundson, and N. Stikov, 2009
%  (c) Board of Trustees, Leland Stanford Junior University    

clear all
close all

T1path = '~/Documents/MATLAB/CEST_EVAL-master-2/T1mapping_060210/';

% Where to find the DICOM data 
loadpath = '/Users/caterina/Desktop/CEST/WISC_backup/Caterina/DATASETS/rats_tumor/week3/MaAm_FOR_14_06_Glio_dynFMISO_Rat4_10012017_1_1_20170110_100236/T1dicoms';

% Where to save the .mat 
savename = [T1path 'TestSingleSlice'];

currpath = pwd;
cd(loadpath)
addpath([T1path '/mfiles']);

d = dicomLoadAllSeries('.');
cd(currpath)

nbrow = size(d(1).imData,1);
nbcol = size(d(1).imData,2);
nbslice = size(d(1).imData, 3);

nbseries = length(d);

data = zeros(nbrow,nbcol,nbslice,nbseries); % Complex data
extra.tVec = zeros(1,nbseries); % One series corresponds to one TI

for k = 1:nbseries
	dataTmp = d(k).imData;
	dataTmp = double(squeeze(dataTmp));	
	for ss = 1:nbslice
		%data(:,:,ss,k) = dataTmp(:,:,3+(ss-1)*4)+i*dataTmp(:,:,4+(ss-1)*4); %SE
		%data(:,:,ss,k) = dataTmp(:,:,1+(ss-1)*4); %GRE
		data(:,:,ss,k) = dataTmp(:,:,ss); 
	end
	extra.tVec(k) = d(k).inversionTime;
end 

% This is where you would correct the phase of certain datapoints. 
% Indeed, the chopping procedure during Prescan can give an additional 180?
% phase on the whole image for certain TIs (make also sure you use manual
% Prescan and don't change gains/center frequency from TI to TI!)
% e.g., data(:,:,3) = -data(:,:,3); 
%    or data(:,:,:,3) = -data(:,:,:,3); if you have more than one slice
% !!! :,:,1 is the first series acquired, not the smallest TI !!!

extra.T1Vec = 1:5000; % this range can be reduced if a priori information is available

data(:,:,2) = -data(:,:,2); 

TI = extra.tVec
save(savename,'data','extra')
