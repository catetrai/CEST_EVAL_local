function [fig, ix] = protocolUI(protocol, sequenceType, ix, fig)
% ** function protocolUI(protocol, sequenceType)
%
% User Interface for loading CEST images. Displays the 'protocol' structure
% as a table and allows the user to select files corresponding to the
% desired 'sequenceType'. 'fig' is the handle to the table figure.
%
% CT 20170318

hasSelectionBeenMade = false;

% figure properties
if nargin<4
    fig = figure('NumberTitle', 'off', 'MenuBar', 'none');
    fig.Position = get(fig, 'Position') .* [1 1 2 1.3];
end
fig.Name = sprintf('Select protocol entry for %s sequence(s)', sequenceType);

% make UI table
t = uitable(fig);
t.Position = t.Position .* [1 1 3.2 1.7];
t.ColumnWidth = {190 100 110 80 70 50 50 50 50 'auto' 'auto' 150 70 'auto'};
uitdata = struct2cell(protocol)';
for i=1:numel(uitdata)
    if isnumeric(uitdata{i}) && ~isscalar(uitdata{i})
        uitdata{i}=num2str(uitdata{i});
    end
end
t.Data = uitdata;
t.ColumnName = fieldnames(protocol)';
t.ColumnFormat = cellstr(repmat('char',[size(uitdata,2) 1]))';
t.RowStriping = 'on';
t.CellSelectionCallback = @getSelectionDirectory;

% make button
btnPos(1:2) = t.Position(1:2) + t.Position(3:4) + [20 -50];
btnPos(3:4) = [100 40];
btn = uicontrol('Style', 'pushbutton', 'String', sprintf('Set %s',sequenceType),...
    'Position', btnPos, 'Callback', @clearFigure);

% wait to clear figure until button push
waitfor(fig, 'Name', '');

% functions
    function getSelectionDirectory(hObject, eventdata)
        ix(1).(sequenceType) = eventdata.Indices(:,1);
        hasSelectionBeenMade = true;
%         assignin('base', ixvar, ind);
    end

    function clearFigure(hObject, eventdata)
%         if evalin('base', sprintf('exist(''%s'',''var'')',ixvar))
        if hasSelectionBeenMade
            clf;
            fig.Name = '';
            pause(0.5);
        else
            error('Please select at least one protocol entry.')
        end
    end
end