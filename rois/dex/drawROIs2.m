function drawROIs(fRefreshRoiTraces)

hIm = findobj(gca, 'Type', 'image');
catchMouseDrag(hIm, @DragStartFun, @CursorMotionFun, @DropFun);

function hP = DragStartFun(hIm, hAx, currentPoint);
    hP = patch(NaN, NaN, 'b', 'FaceColor', 'none', 'EdgeColor', [1 1 1]*0, 'LineWidth', 2, 'Tag', 'ROIPatch');
end
    
function CursorMotionFun(hIm, hAx, mouseTrace, hP);
    set(hP, 'XData', [mouseTrace(:, 1); mouseTrace(1, 1)], 'YData', [mouseTrace(:,2); mouseTrace(1,2)])
end
    
function DropFun(hIm, hAx, mouseTrace, hP);
    setupPatchBehavior2(hP, fRefreshRoiTraces);
end
end