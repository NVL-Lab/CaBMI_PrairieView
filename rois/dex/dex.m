function fAccessvars = dex(matin,framerate,duration);

%%
dimsizes = size(matin)
%roiPatchHandles = cell([1 1 1 1 dimsizes(5:end)]);
roiPatchHandles = cell([1 1 1 1 1 dimsizes(6:end)]);

h = buildUI(dimsizes);

set(h.Cf, 'Name', ['', inputname(1)])
set([h.Bgx h.Bgy], 'SelectionChangeFcn', @callback)
set(h.Play, 'Callback', @(hObject, varargin) play(hObject, h))
set(h.Sl, 'Callback', @callback)

callback(NaN);
refreshImage
axes(h.Ax);
drawROIs2(@plotRoiTraces)

%keyboard
%uiwait

fAccessvars = @accessvars;

    function out = accessvars(varname)
        out = eval(varname);
        disp(['accessing ', varname])
    end

    function refreshImage(varargin)
        %axis(h.Ax, 'equal', 'tight')
        axis(h.Ax,'tight');
        axes(h.Ax), %autoscale
    end

    function roiTraces=plotRoiTraces;
        if ~(size(matin, 3) > 1), return, end
        
        timedim = 3;
        subsstr = getSubsStruct(h);
        subsstr.subs{timedim} = ':';
        mat = subsref(matin, subsstr);
        
        hP = findobj(h.Ax, 'Tag', 'ROIPatch');
        subsstr.subs([size(roiPatchHandles) ones(1, 10)]==1) = {':'};
        roiPatchHandles = subsasgn(roiPatchHandles, subsstr, {hP});
        
        roiTraces = roitrace(mat, extractROIs2(ancestor(h.Ax, 'figure')));
        roiTraces = permute(roiTraces(:,:,1:end,:,:,:,:,:), [3 1 2 4:10]);
        
        if ~ishandle(h.RoiPlotAx)
            plotfig = figure; plotfigpos = get(plotfig, 'Position');
            imfig = ancestor(h.Ax, 'figure');
            iptwindowalign(imfig,'top',plotfig,'top');
            iptwindowalign(imfig,'right',plotfig,'left');
            h.RoiPlotAx = axes;

%             %This part was added by Emre on 07.06.11
%             framerate=0.197 ;
%             duration=0.197 * 10000;
            set(gca,'XTick',(0:1/framerate:duration/framerate))
            set(gca,'XTickLabel',(0:1:duration));
            xlabel('time (seconds)','fontsize',8,'fontweight','b');
            ylabel('DF/F or F','fontsize',8,'fontweight','b');
            box off; set(gca,'TickDir','out','fontsize',6);
%             save('ROI_' ,roiTraces);
        end
        
        hP = findobj(h.Ax, 'Tag', 'ROIPatch');
        for i = 1:size(hP,1)
            if ~ishandle(getUD(hP(i), 'hPlotTrace'))
                hold(h.RoiPlotAx, 'on')
                hPlotTrace = plot(h.RoiPlotAx, roiTraces(:, i));
                hold(h.RoiPlotAx, 'off')
                setUD(hP(i), 'hPlotTrace', hPlotTrace)
                
                co = get(h.RoiPlotAx, 'ColorOrder');
                currentcolor = 0.8*co(mod(find(sort(hP) == hP(i))-1, size(co, 1))+1, :);
                
                if get(hP(i), 'EdgeColor') == [0 0 0]
                    set(hP(i), 'EdgeColor', currentcolor)
                end
                set(hPlotTrace, 'Color', get(hP(i), 'EdgeColor'))
            else
                set(getUD(hP(i), 'hPlotTrace'), 'YData', roiTraces(:, i))
            end
        end
    end

    function callback(hObject, varargin)
        switch(hObject)
            case num2cell(h.Sl) %any slider handle
                val = round(get(hObject, 'Value'));
                set(hObject, 'Value', val, 'TooltipString', num2str(val));
                set(h.CurrSliceLabel(h.Sl == hObject), 'String', num2str(val));
                if find(h.Sl == hObject) ~= 3
                    hP = findobj(h.Ax, 'Tag', 'ROIPatch');
                        for i = 1:size(hP,1)
                            if ishandle(getUD(hP(i), 'hPlotTrace'))
                                delete(getUD(hP(i), 'hPlotTrace'))
                            end
                            set(hP(i), 'Visible', 'off', 'HandleVisibility', 'off') %delete(hP(i))
                        end

                    subsstr = getSubsStruct(h); subsstr.type = '{}'; subsstr.subs([size(roiPatchHandles) ones(1, 10)]==1) = {':'};
                    set(subsref(roiPatchHandles, subsstr)','Visible', 'on', 'HandleVisibility', 'on')
                    plotRoiTraces
                end
            case {h.Bgx, h.Bgy}
                displayimage
                refreshImage
                plotRoiTraces
                
        end
        displayimage
        
    end
    
    function displayimage
        im = squeeze(subsref(matin, getSubsStruct(h)));

        xydim = getXyDim(h);
        if xydim(2) < xydim(1), im = im'; end
        
        %im = imfilter(im, fspecial('gaussian', [3 3], 0.5));
        
        set(h.Im, 'CData', im);
        drawnow
        %autoscale;
    end

    function ImButtonDown(hObject, varargin)
        cp = get(h.Ax, 'CurrentPoint');
        cp = round(cp(1, 1:2));
        disp(cp(2))
        if ~ishandle(h.SlicePlot)
            h.SlicePlot = figure;
        end
        xydim = getXyDim(h);
        sb = getSubsStruct(h, xydim(2));
        sb.subs{xydim(1)} = cp(2);
        y = subsref(matin, sb);
        y = double(squeeze(y));
        %keyboard
        hPlotAx = findobj(h.SlicePlot, 'Type', 'axes');
        plot(hPlotAx, y(1:end))
        
    end

    function play(hObject, h)
        fps = 10;
        sliderhandle = h.Sl(get(hObject, 'UserData'));
        maxframes = get(sliderhandle, 'Max');
        tic
        while (get(hObject, 'Value') == 1)
            val = round(get(sliderhandle, 'Value')+1);
            if val > maxframes, set(sliderhandle, 'Value', 1), displayimage, break, end
            set(sliderhandle, 'Value', val)
            displayimage
            pause(1/fps - toc), tic
        end
        set(hObject, 'Value', 0)
    end

end

        



function h = buildUI(dimsizes, fCallback, fPlay, fRefreshImage)
    h.RoiPlotAx = NaN;
    h.SlicePlot = NaN;
    h.Imfig = figure('Position', [0 0 512 512], 'Visible', 'off');
    h.Ax = axes('Parent', h.Imfig); set(h.Ax, 'Position', [0 0 1 1])
    movegui(h.Imfig, 'northwest'); imfigpos = get(h.Imfig, 'Position');
    set(h.Imfig, 'Visible', 'on')
    h.Cf = figure;
    h.Im = imagesc(NaN, 'Parent', h.Ax, 'EraseMode', 'normal');
    
    
    %fpos = get(h.Cf, 'Position'); 
    fpos = [imfigpos(1) imfigpos(2)-80-length(dimsizes)*20 400 length(dimsizes)*20+5];
    set(h.Cf, 'Position', fpos, 'MenuBar', 'none', 'Name', 'dex', 'Resize', 'off', 'NumberTitle', 'off');
    %iptwindowalign(h.Imfig, 'bottom', h.Cf, 'top'), iptwindowalign(h.Imfig, 'left', h.Cf, 'left')
    %set(h.Cf, 'WindowButtonUpFcn', @(varargin) set(h.Cf, 'WindowButtonMotionFcn', ''))
    set(h.Cf, 'BusyAction', 'cancel', 'WindowButtonMotionFcn', @cursorMotion)
    
    h.Tb = uitoolbar(h.Cf);
    h.tbRois = uitoggletool('CData', rand(20,20,3));
    h.tbAutosc = uipushtool('CData', rand(20,20,3), 'ClickedCallback', @(varargin) autoscale);
    
    h.Bgx = uibuttongroup('BorderType', 'none', 'HitTest', 'off');
    h.Bgy = uibuttongroup('BorderType', 'none', 'HitTest', 'off');
    for iDim = 1:length(dimsizes)        
        ds = dimsizes(iDim);
        h.Rbx(iDim) = uicontrol(h.Cf, 'Style', 'Radiobutton', 'Position', [5 fpos(4)-20*(iDim) 15 15], 'Parent', h.Bgx, 'UserData', iDim);
        h.Rby(iDim) = uicontrol(h.Cf, 'Style', 'Radiobutton', 'Position', [20 fpos(4)-20*(iDim) 15 15], 'Parent', h.Bgy, 'UserData', iDim);
        h.Sl(iDim) = uicontrol(h.Cf, 'Style', 'Slider', ...
            'SliderStep', [1/(max([ds 2])-1) 1/(max([ds 2])-1)], 'Value', 1, 'Min', 1, 'Max', dimsizes(iDim), ...
            'Position', [130 fpos(4)-20*(iDim)+1 230 15], 'HitTest', 'on');
        h.DimLabel(iDim) = uicontrol(h.Cf, 'Style', 'text', 'Position', [35 fpos(4)-20*(iDim) 40 15], ...
            'String', ['dim: ', num2str(iDim)], 'HorizontalAlignment', 'left');
        h.CurrSliceLabel(iDim) = uicontrol(h.Cf, 'Style', 'text', 'Position', [70 fpos(4)-20*(iDim) 30 15], ...
            'String', ['1'], 'HorizontalAlignment', 'right');
        h.SliceNumLabel(iDim) = uicontrol(h.Cf, 'Style', 'text', 'Position', [100 fpos(4)-20*(iDim) 30 15], ...
            'String', ['/', num2str(dimsizes(iDim))], 'HorizontalAlignment', 'left');
        h.Play(iDim) = uicontrol(h.Cf, 'Style', 'togglebutton', 'Position', [365 fpos(4)-20*(iDim)+1 30 15], 'String', 'Play', 'UserData', iDim);
        %set(hLabel(iDim), 'BackgroundColor', get(h.Cf, 'Color'))
    end 
    set(h.Bgx,'SelectedObject',h.Rbx(1))
    set(h.Bgy,'SelectedObject',h.Rby(2))

    
    function cursorMotion(hCf, varargin)
        if isempty(gco) | ~any(h.Sl == gco), return, end
        persistent oldvalue, if isempty(oldvalue), oldvalue = 0; end        
        val = round(get(gco, 'Value'));
        if val ~= oldvalue
            slcallback = get(gco, 'Callback');
            slcallback(gco);
            oldvalue = val;
        end           
    end

end

function subsStruct = getSubsStruct(h, xydim)
    if nargin < 2, xydim = getXyDim(h); end
    subsStruct.type='()';
    subsStruct.subs = num2cell(getSelectedDims(h));
    subsStruct.subs(xydim) = {':'};
end

function selectedDims = getSelectedDims(h)
    selectedDims = (round(cell2mat(get(h.Sl, 'Value'))));
end

function xydim = getXyDim(h)
    xydim(1) = get(get(h.Bgx,'SelectedObject'), 'UserData');
    xydim(2) = get(get(h.Bgy,'SelectedObject'), 'UserData');
end




% %%
% import java.awt.*
% import javax.swing.*
% import com.mathworks.mwt.*
% 
% dimsize = [2 4 2];
% 
% frame = JFrame('title')
% frame.addWindowListener(window.MWWindowActivater(frame));
% pane = frame.getContentPane();
% pane.setLayout(GridLayout(3,2));
% 
% for iDim = 1:length(dimsize)
%     label = JLabel(num2str(dimsize(iDim)))
%     slider(iDim) = JSlider;
%     slider(iDim).setPaintTicks(true);
%     slider(iDim).setSnapToTicks(true);
%     slider(iDim).setPaintLabels(true);
%     
%     slider(iDim).setMinimum(1);
%     slider(iDim).setMaximum(dimsize(iDim));
%     slider(iDim).setMinorTickSpacing(1);
%     slider(iDim).setMajorTickSpacing(min([5 dimsize(iDim)-1]));
%     pane.add(label)
%     pane.add(slider(iDim))
% end
% 
% frame.pack();
% frame.setVisible(true);