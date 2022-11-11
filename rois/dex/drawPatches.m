function drawPatches(hAx, patchCoords)

if ~iscell(patchCoords),
    if isnan(patchCoords), return, end
else
    if isempty(patchCoords{1}), return, end
end

for i = fliplr(1:size(patchCoords, 1))
    hP = patch(patchCoords{i}(:, 1), patchCoords{i}(:, 2), ...
        'b', 'FaceColor', 'none', 'EdgeColor', [1 1 1]*0, 'LineWidth', 2, 'Tag', 'ROIPatch', 'Parent', hAx);
end