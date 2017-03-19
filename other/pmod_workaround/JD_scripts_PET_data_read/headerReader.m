%HEADERREADER opens a Siemens Inveon (R) header file and returns all
%       variables in a structure.
%
%   USAGE: 
%       headerReader;
%           Will start with a dialog box to select a header file.
%       headerReader(fileName);
%           Will process the headerfile given in 'fileName'
%
% J.A. Disselhorst, 2009-2016
% Univeristy of Twente, Enschede (NL)
% Radboud University Medical Center, Nijmegen (NL)
% Werner Siemens Imaging Center, Tuebingen (DE)
%
% Version 2016.01.10
%      Last version: bug in frames fixed.
%
% Disclaimer:
% THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT WARRANTY OF ANY
% KIND, EITHER EXPRESSED OR IMPLIED AND IS TO BE USED AT YOUR OWN RISK 

function [data,fileName] = headerReader(fileName)

if nargin~=1
    [fileName,PathName] = uigetfile({'*.hdr','Header files (*.hdr)'},'Select header file');
    if fileName
        fileName = fullfile(PathName,fileName);
    else
        fprintf('no file selected\n')
        return
    end
end

try
    fid = fopen(fileName);
    dataout = textscan(fid,'%s','delimiter','\n','commentstyle','#');
    dataout = dataout{:};
    fclose(fid);
catch %#ok<CTCH>
    error(['File: ''' fileName ''' does not exist or cannot be opened.']);
end
nLines = size(dataout,1);
secStr = 'General'; sectionCount = -1;
data = struct(secStr,[]);

for ii=1:nLines
    line = cell2mat(dataout(ii));
    ind = find( ~isspace(line) );       % indices of the non-space characters in the str    
    if isempty(ind)
        line = [];        
    else
        line = line( ind(1):ind(end) );
    end
    
    value = [];
    key = [];
    if isempty(line)                            % empty line
        status = 0; 
    else                                        % key found, perhaps also value.
        pos = find(isspace(line),1,'first');
        if ~isempty(pos)                        % key-value pair found
            status = 2;
            key = lower(line(1:pos-1));
            value = line(pos+1:end);
            
            if all(~strcmp(key,{'singles','projection','ct_projection_center_offset','ct_projection_horizontal_bed_offset','projection_duration','spect_pinhole_params'}))
                if isempty(regexp(value,'[^\d\s.-+e]','once'))
                    temp = str2double(regexp(value, '-*\w*[.]*\w*[e+-]*\w*', 'match'));
                    if ~any(isnan(temp)) && length(regexp(value,'[.]'))<2     %% If result is a valid number. i.e. max one dot. etc.
                        value = temp;
                    end
                end
            end
            
            if isempty(key)                     % empty keys are not allowed
                status = 0;
                key = [];
                value = [];
            end
        else                                    % empty value
            status = 3;
            key = line;
            if strcmp(key,'end_of_header')
                if ii<nLines
                    sectionCount = sectionCount+1;
                    value = sprintf('frame_%u',sectionCount);
                    status = 1;
                else
                    return
                end
            end
        end
    end

    % Probably not necessary:
    try
        key = regexprep(key,'\W','_'); %replace illegal characters in keynames with underscore
    end
    if length(key)>1
        key = [regexprep(key(1),'[\d_]','X$0') key(2:end)]; %find underscores and numbers at first position and add an 'X'
    elseif length(key)==1
        key = regexprep(key(1),'[\d_]','X$0');  %find underscores and numbers at first position and add an 'X'
    end

        

    if status == 1 % new section
        secStr = value;
        data.(secStr)=[];
    elseif status == 2 % new key-value pair
        if ~isfield( eval(['data.' secStr]), key )
            data.(secStr).(key) = value;
        else
            temp = getfield(data,secStr,key);
            if isnumeric(temp) && isnumeric(value)   %only numbers
                temp = [temp; value];
            elseif ~iscell(temp)                    %strings and stuff
                temp = cellstr(temp);
                temp(end+1,1) = {value};
            else
                temp(end+1,1) = {value};
            end
            data.(secStr).(key) = temp;
        end
    elseif status == 3 % empty key
        data.(secStr).(key) = ' ';
    end
end
if isempty(data.(secStr))
    data = rmfield(data,secStr);
end


