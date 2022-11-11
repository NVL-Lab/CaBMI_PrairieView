function setupPatchBehavior(hPin, fRefreshRoiTraces)

for i = 1:size(hPin, 1)
    setUD(hPin(i), 'hPlotTrace', NaN);

    if size(get(hPin(i), 'XData'), 1) < 5,
        delete(hPin(i))
        return
    end

    set(hPin(i), 'UICOntextMenu', makeCmRoi(hPin(i)));
    set(hPin(i), 'ButtonDownFcn', @patchClick);

    CurrentPoint = getCurrentPoint(get(hPin(i), 'Parent'));
    catchMouseDrag(hPin(i), @patchClick, @patchMotionFun, @patchDropFun, 'extend');
    refreshROIs(hPin(i))
end
    fRefreshRoiTraces()
    
    function patchMotionFun(hP, hAx, mouseTrace, custom);
        initxdata = get(hP, 'XData');
        initydata = get(hP, 'YData');
        dist = mouseTrace(end,:) - mouseTrace(end-1,:);
        set(hP, 'XData', initxdata + dist(1))
        set(hP, 'YData', initydata + dist(2))
        refreshROIs(hP)

        fRefreshRoiTraces()
    end


    function patchDropFun(hP, hAx, mouseTrace, custom);
        if ishandle(getUD(hP, 'hPlotTrace'))
            set(getUD(hP, 'hPlotTrace'), 'LineWidth', 0.5)
        end

        refreshROIs(hP)
        fRefreshRoiTraces()
    end

    function refreshROIs(hP)
        hAx = ancestor(hP, 'axes');
        hIm = findobj(hAx, 'Type', 'Image');
        s = size(get(hIm, 'CData'));
        xdata = get(hP, 'XData');
        ydata = get(hP, 'YData');
        binroi = poly2mask(xdata, ydata, s(1), s(2));
        %set(hP, 'UserData', binroi);
        setUD(hP, 'binroi', binroi);
    end

    function hCmRoi = makeCmRoi(hP)
        hCmRoi = uicontextmenu('Parent', ancestor(hP, 'figure'));
        uimenu(hCmRoi, 'Label', 'Delete', 'Callback', {@CmRoiCallback});
        uimenu(hCmRoi, 'Label', 'Copy', 'Callback', {@CmRoiCallback});
        uimenu(hCmRoi, 'Label', 'Delete All', 'Callback', {@CmRoiCallback}, 'Separator', 'on');
        uimenu(hCmRoi, 'Label', 'Highlight all', 'Callback', {@CmRoiCallback});
    end

    function CmRoiCallback(hObj, varargin)
    hF = ancestor(hObj, 'figure');
    hCmRoi = get(hObj, 'Parent');
    hP = findobj(hF, 'UICOntextMenu', hCmRoi);
    hPall = findobj(hF, 'Tag', 'ROIPatch');

    switch get(hObj, 'Label')
        case 'Highlight all'
            switch(get(hObj, 'checked'))
                case 'off'
                    set(hObj, 'checked', 'on')
                    set(hPall, 'FaceColor', [1 1 1]*0)
                case 'on'
                    set(hObj, 'checked', 'off')
                    set(hPall, 'FaceColor', 'none')
            end
        case 'Delete'
            if ishandle(getUD(hP, 'hPlotTrace'))
                delete(getUD(hP, 'hPlotTrace'))
            end

            delete(hP)
        case 'Copy'

            hAx = get(hP, 'Parent');
            hCopy = copyobj(hP, hAx);
            disp(2334)
            setUD(hCopy, 'hPlotTrace', NaN);
            set(gcf, 'SelectionType', 'extend');
            fButtonDown = get(hCopy, 'ButtonDownFcn');
            fButtonDown(hCopy)
            set(hCopy, 'UICOntextMenu', makeCmRoi(hP));
            refreshROIs(hCopy)
        case 'Delete All'
            hAx = get(hP, 'Parent');
            for i = 1:size(hPall,1)
                if ishandle(getUD(hPall(i), 'hPlotTrace'))
                    delete(getUD(hPall(i), 'hPlotTrace'))
                end            
                delete(hPall(i))
            end
    end

    fRefreshRoiTraces()
    end
end


function CurrentPoint = patchClick(hP, hAx, varargin)
    hP
    if ishandle(getUD(hP, 'hPlotTrace'))
        set(getUD(hP, 'hPlotTrace'), 'LineWidth', 2)
        initWindowButtonUpFcn = get(ancestor(hP, 'figure'), 'WindowButtonUpFcn');
        set(ancestor(hP, 'figure'), 'WindowButtonUpFcn', @(hF, varargin) WindowButtonUpFcn(hF, hP))
    end
    function WindowButtonUpFcn(hF, hP)
        %if ~isempty(initWindowButtonUpFcn), initWindowButtonUpFcn(hF, varargin), end
        set(getUD(hP, 'hPlotTrace'), 'LineWidth', 0.5)
        %set(ancestor(hP, 'figure'), 'WindowButtonUpFcn', initWindowButtonUpFcn)
    end
    
    CurrentPoint = getCurrentPoint(get(hP, 'Parent'));
end