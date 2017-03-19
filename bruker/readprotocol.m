function P = readprotocol(directory)
% ** function P = readprotocol(directory)
%
% Given patient directory, creates structure 'P' with meta-data for each
% sequence in the study protocol.
% 'P' has the following fields:
%    'ScanName': name of the sequence set by the user at the scanner
%    'ProtocolName'
%    'Method'
%    'NRepetitions': number of instances/acquisitions within the sequence
%    'SliceOffset': in mm
%    'B1': saturation pulse strength (in microTesla)
%    'TI': inversion times, if Selective Inversion Recovery sequence
%    'TR'
%    'TE'
%    'MatrixSize'
%    'FOV': in mm
%    'SpatResol': in mm
%    'SliceThick': in mm
%    'Path': full path of the sequence directory
%
% CT 20170312

P = dir(directory);
nonNumericMask = cellfun('isempty', regexp({P.name},'^\d+$'));
P(nonNumericMask) = [];  % sequence folder names must be numbers
folderNum = str2double({P.name}');

[P.Path_temp] = P.name;
P = rmfield(P, {'name', 'date', 'bytes', 'isdir', 'datenum'});

h = waitbar(0, 'Loading study protocol...');
for i=1:length(P)
    P(i).Path_temp = fullfile(directory, P(i).Path_temp);
    
    method_path = fullfile(P(i).Path_temp, 'method');
    acqp_path = fullfile(P(i).Path_temp, 'acqp');
    try
        P(i).ScanName = mineMetaDataFile(acqp_path, 'ACQ_scan_name');
        P(i).ScanName = P(i).ScanName(2:end-1);
    catch
        P(i).ScanName = '';
    end
    try
        P(i).ProtocolName = mineMetaDataFile(acqp_path, 'ACQ_protocol_name');
        P(i).ProtocolName = P(i).ProtocolName(2:end-1);
    catch
        P(i).ProtocolName = '';
    end
    try
        P(i).Method = mineMetaDataFile(acqp_path, 'ACQ_method');
        P(i).Method = P(i).Method(2:end-1);
    catch
        P(i).Method = '';
    end
    try
        P(i).NRepetitions = mineMetaDataFile(method_path, 'PVM_NRepetitions');
    catch
        P(i).NRepetitions = NaN;
    end
    try
        P(i).SliceOffset = mineMetaDataFile(method_path, 'PVM_SliceOffset');
    catch
        P(i).SliceOffset = NaN;
    end
    try
        if strcmp(mineMetaDataFile(method_path, 'PVM_MagTransOnOff'), 'On')
            try
                P(i).B1 = mineMetaDataFile(method_path, 'PVM_MagTransPower');
            catch
                P(i).B1 = NaN;
            end
        else
            P(i).B1 = NaN;
        end
    catch
        P(i).B1 = NaN;
    end
    try
        if strcmp(mineMetaDataFile(method_path, 'PVM_SelIrOnOff'), 'On')
            try
                P(i).TI = mineMetaDataFile(method_path, 'PVM_SelIrInvTime');
            catch
                P(i).TI = NaN;
            end
        else
            P(i).TI = NaN;
        end
    catch
        P(i).TI = NaN;
    end
    try
        P(i).TR = mineMetaDataFile(method_path, 'PVM_RepetitionTime');
    catch
        P(i).TR = NaN;
    end
    try
        P(i).TE = mineMetaDataFile(method_path, 'PVM_EchoTime');
    catch
        P(i).TE = NaN;
    end
    try
        P(i).MatrixSize = mineMetaDataFile(method_path, 'PVM_Matrix');
    catch
        P(i).MatrixSize = NaN;
    end
    try
        P(i).FOV = mineMetaDataFile(method_path, 'PVM_Fov');
    catch
        P(i).FOV = NaN;
    end
    try
        P(i).SpatResol = mineMetaDataFile(method_path, 'PVM_SpatResol');
    catch
        P(i).SpatResol = NaN;
    end
    try
        P(i).SliceThick = mineMetaDataFile(method_path, 'PVM_SliceThick');
    catch
        P(i).SliceThick = NaN;
    end
    P(i).Path = P(i).Path_temp;
    
    waitbar(i/length(P), h);
end
P = rmfield(P, {'Path_temp'});

% sort entries by increasing folder number (to preserve sequence order)
[~,ixsort] = sort(folderNum);
P = P(ixsort);

close(h);