function rois = extractROIs(hF)

rois = [];

if nargin < 1
    hF = gcf;
end

hP = findobj(hF, 'Tag', 'ROIPatch');
for i = 1:size(hP, 1)
    rois(:,:,i) = getUD(hP(i), 'binroi');
end
%if ~iscell(rois), rois = {rois}; end
%rois = cat(3, rois{:});