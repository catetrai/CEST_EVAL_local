%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD_CEST_BRUKER %%%%%%%%%%%%%%%%%%%%%%%%%%%%

USEGUI = true;  % use User Interface for selecting files? (recommended)
if ~exist('READIMAGES', 'var')
    READIMAGES = true;    % read in the images as default
end

%% Get patient directory & study protocol
if evalin('base', '~exist(''directory'', ''var'');');
    directory = uigetdir('', 'Select study directory');
    protocol = readprotocol(directory);
end

%% Select scans from protocol
% the 'ix' structure will contain indices of selected protocol entries
if ~exist('ix', 'var')
    ix = struct([]);
end

if USEGUI
    [uifig, ix] = protocolUI(protocol, sequenceType, ix);
    if any(strcmp(sequenceType, {'Mz', 'WASABI'}))
        [~, ix] = protocolUI(protocol, 'M0', ix, uifig);
    end
    close(uifig);
else
    inputmsg = sprintf('\nEnter protocol entry number(s) for %s sequence(s): ', sequenceType);
    ix(1).(sequenceType) = input(inputmsg); clear inputmsg
    if any(strcmp(sequenceType, {'Mz', 'WASABI'}))
        x(1).M0 = input('\nEnter protocol entry number for M0 sequence: ');
    end
end

%% Read images as arrays of doubles & write scan parameters into 'P' structure
if READIMAGES
    switch sequenceType
        case {'Mz', 'WASABI'}
            Mz_stack = load_Mz(protocol(ix.(sequenceType)).Path);   % dimensions: Mz_stack(x,y,1,w)
            M0_stack = load_M0(protocol(ix.M0).Path);               % dimensions: M0_stack(x,y)
            
            % Frequency offets (deltaomega in ppm) are stored in 'P.SEQ.w'
            P = wipread_modified(protocol(ix.(sequenceType)).Path, protocol(ix.M0).Path);
        case 'T1mapping'
            T1_stack = load_T1(protocol(ix.(sequenceType)).Path);   % dimensions: T1_stack(x,y,1,nTI)
    end
end