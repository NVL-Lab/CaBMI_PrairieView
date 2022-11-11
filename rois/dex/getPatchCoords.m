function patchCoords = getPatchCoords(hAx)

if nargin < 1
    hAx = gca;
end


hP = findobj(hAx, 'Tag', 'ROIPatch');

xdata = get(hP, 'XData');
ydata = get(hP, 'YData');

if ~iscell(xdata), xdata = {xdata}; ydata = {ydata}; end

for i = 1:size(xdata, 1)
    patchCoords{i, :} = [xdata{i} ydata{i}];
end