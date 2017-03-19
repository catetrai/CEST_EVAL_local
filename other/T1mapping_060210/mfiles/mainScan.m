% mainScan.m
%
% Loads .mat file from the directory 'data' in the T1path, 
% performs T1 mapping, and displays the results 
%
% written by J. Barral, M. Etezadi-Amoli, E. Gudmundson, and N. Stikov, 2009
%  (c) Board of Trustees, Leland Stanford Junior University

clear all
close all

T1path = '~/Documents/MATLAB/CEST_EVAL-master-2/T1mapping_060210/';

%% Where to find the data
loadpath = '/Users/caterina/Desktop/CEST/WISC_backup/Caterina/DATASETS/rats_tumor/week3/MaAm_FOR_14_06_Glio_dynFMISO_Rat4_10012017_1_1_20170110_100236/T1dicoms';

datasetnb = 1;

switch (datasetnb)
	case 1
		filename = 'TestSingleSlice'; % complex fit
		method = 'RD-NLS'
	case 2
		filename = 'TestSingleSlice'; % magnitude fit
		method = 'RD-NLS-PR'
end


%% Where to save the data
savepath = [T1path 'fitdata/'] 

loadStr = [T1path filename]
saveStr = [savepath 'T1Fit' method '_' filename]

%% Perform fit
T1ScanExperiment(loadStr, saveStr, method);

%% Display results
T1FitDisplayScan(loadStr, saveStr);
