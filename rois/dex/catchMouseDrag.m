function [mouseTrace] = catchMouseDrag(hObj, fDragStart, fDragMotion, fDrop, SelectionType)
%CATCHMOUSEDRAG function to catch mouse movement
%
%   [mouseTrace] = catchMouseDrag(hObj, fDragStart, fDragMotion, fDrop, SelectionType)
%   
%   Sets up mouse movement listeners to catch cursor movement within a
%   MATLAB figure. 
%
%   <hObj>          Handle of the object to draw on
%   <fDragStart>    Handle to the function called when the mouse button is
%                   pressed. Syntax: fDragStart(hIm, hAx, currentPoint)
%   <fDragMotion>   Handle to the function called when the mouse cursor is
%                   moved while the button is still pressed. Syntax: 
%                   fDragMotion(hIm, hAx, mouseTrace, hP)
%   <fDrop>         Handle to the function called when the mouse button is
%                   released. Syntax: fDrop(hObj, hAx, mouseTrace, customData);
%   <SelectionType> 'normal' (default), 'extend', 'alt' or 'open'. See
%                   figure SelectionType property for more information.
%
%   Example 1:
%     figure;
%     hIm = image;
%     mousetrace = catchMouseDrag(gca) %now draw a scribble
%
%   Example 2:
%     %file drawROIs.m start:
%     function drawROIs
%     hAx = gca;
%     hIm = findobj(hAx, 'Type', 'image');
%     catchMouseDrag(hIm, @DragStartFun, @CursorMotionFun);
%     function hP = DragStartFun(hIm, hAx, currentPoint);
%         hP = patch(NaN, NaN, 'b', 'FaceColor', 'none', 'EdgeColor', [1 1 1]*0, 'LineWidth', 2.5, 'Tag', 'ROIPatch', 'Parent', hAx);  
%     function CursorMotionFun(hIm, hAx, mouseTrace, hP);
%         set(hP, 'XData', [mouseTrace(:, 1); mouseTrace(1, 1)], 'YData', [mouseTrace(:,2); mouseTrace(1,2)])
%     %file drawROIs.m end

%   2005: created, BJ
%   060922: updated to remove conflicts between successive calls, BJ 


    if nargin < 5
        SelectionType = 'normal';
    end, if nargin < 4
        fDrop = @(varargin) NaN;
    end, if nargin < 3
        fDragMotion = @(varargin) NaN;
    end, if nargin < 2
        fDragStart = @(varargin) NaN;
    end
    
    initButtonDownFcn = get(hObj, 'ButtonDownFcn');
    set(hObj, 'ButtonDownFcn', ...
        @(hObj, ignore) objButtonDownFcn(hObj, fDragStart, fDragMotion, fDrop, SelectionType, initButtonDownFcn));
    if nargout > 0, uiwait, end
    
    function objButtonDownFcn(hObj, fDragStart, fDragMotion, fDrop, SelectionType, initButtonDownFcn)
        
        if ~isempty(initButtonDownFcn), initButtonDownFcn(hObj);, end
        hF = ancestor(hObj, 'figure');
        hAx = ancestor(hObj, 'axes');
        initWindowButtonMotionFcn = ''; %get(hF, 'WindowButtonMotionFcn');
        initWindowButtonUpFcn = ''; %get(hF, 'WindowButtonUpFcn');
        if ~strcmp(get(hF, 'SelectionType'), SelectionType), return, end
        mouseTrace = getCurrentPoint(hAx);
        customData = fDragStart(hObj, hAx, mouseTrace);
        set(hF, 'WindowButtonMotionFcn', {@figWindowButtonMotionFcn, customData, hObj}, ...
            'WindowButtonUpFcn', @figWindowButtonUpFcn)

        function figWindowButtonUpFcn(hF, varargin)
            set(hF, 'WindowButtonMotionFcn', initWindowButtonMotionFcn, ...
                'WindowButtonUpFcn', initWindowButtonUpFcn);
            fDrop(hObj, hAx, mouseTrace, customData);
            uiresume
        end

        function figWindowButtonMotionFcn(hF, varargin, customData, hObj)
            currentPoint = getCurrentPoint(hAx);
            mouseTrace(end+1,:) = currentPoint;
            fDragMotion(hObj, hAx, mouseTrace, customData);
        end
    end
end