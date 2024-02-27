%% cb_dataTipCursor
%
% Callback function to manage the behavior of dataCursor mode.
% By default, when the cursor in placed over a node of a plotted patch,
% only the coordinates are provided as the tip.
% This function is an adapted version of the default built-in callback
% function to show also the result value that originates the patch plot.
%
% Refs.:
% https://www.mathworks.com/matlabcentral/answers/303793-data-cursor-showing-color-value
% https://www.mathworks.com/help/matlab/ref/matlab.graphics.shape.internal.datacursormanager.html
% https://www.mathworks.com/help/matlab/creating_plots/create-custom-data-tips.html
%
function txt = cb_dataTipCursor(~,info)
pos = info.Position;
dim = length(pos);

if (dim == 2)
    % Coordinates
    x = pos(1);
    y = pos(2);
    
    % nodal value
    xIdx = find(info.Target.XData == x);
    yIdx = find(info.Target.YData == y);
    idx  = intersect(xIdx,yIdx);
    val  = info.Target.CData(idx(1));
    
    % Displaying text
    txt = {['X',formatValue(x,info)],...
           ['Y',formatValue(y,info)],...
           ['Value',formatValue(val,info)]};
    
elseif (dim == 3)
    % Coordinates
    x = pos(1);
    y = pos(2);
    z = pos(3);
    
    % nodal value
    xIdx = find(info.Target.XData == x);
    yIdx = find(info.Target.YData == y);
    zIdx = find(info.Target.ZData == z);
    idx  = intersect(intersect(xIdx,yIdx),zIdx);
    val  = info.Target.CData(idx(1)); 
    
    % Displaying text
    txt = {['X',formatValue(x,info)],...
           ['Y',formatValue(y,info)],...
           ['Z',formatValue(z,info)],...
           ['Value',formatValue(val,info)]};
end

function formattedValue = formatValue(value,info)
if strcmpi(info.Interpreter,'tex')
    valueFormat = ' \color[rgb]{0 0.6 1}\bf';
    removeValueFormat = '\color[rgb]{.25 .25 .25}\rm';
else
    valueFormat = ': ';
    removeValueFormat = '';
end
formattedValue = [valueFormat num2str(value,4) removeValueFormat];
