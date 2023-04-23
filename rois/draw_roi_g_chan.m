function [roi_data] = draw_roi_g_chan(plot_images, roi_data)
%TODO: DATA FOR TWO CHANNELS:
%Allows user to draw shapes onto an image
%roi_data fields: 
% data for pooling over rois
% data for each channel separate
%
% current: 
% --
% num_rois
% roi_bin_cell
% roi_mask
% roi_mask_bin
% im_roi
%
% next: 

% roi_data initialization:
%
% roi_data.im_bg          = im_bg; 
% roi_data.im_roi         = im_bg; 
% roi_data.im_roi_rg      = im_bg; 
% roi_data.num_rows       = size(im_bg,1);
% roi_data.num_cols       = size(im_bg,2);
% roi_data.num_rois       = 0;  
% roi_data.roi_bin_cell   = {};
% roi_data.roi_mask       = zeros(roi_data.num_rows , roi_data.num_cols); 
% roi_data.roi_mask_bin   = zeros(roi_data.num_rows , roi_data.num_cols); 
% roi_data.chan_logical   = []; %num_chan x num_roi
% 
% roi_data.chan = repmat(struct(...
%     'label', '', ...
%     'num_rois', '', ...
%     'idxs', '', ...
%     'im_roi',       zeros(roi_data.num_rows, roi_data.num_cols), ...
%     'roi_mask',     zeros(roi_data.num_rows, roi_data.num_cols), ...
%     'roi_mask_bin', zeros(roi_data.num_rows, roi_data.num_cols)), [2 1]); 



%%
screen_size = get(0,'ScreenSize');
% Show the red and green channels:
close all;
for plot_i=1:length(plot_images)
    im_plot = plot_images(plot_i).im; 
    im_title = plot_images(plot_i).label; 
    h = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);
    imagesc(im_plot); axis square; colormap bone; title(im_title);    
end

disp('Adding ROIs to image!'); 
roi_complete_bool = 0; 
while(~roi_complete_bool)
    %Close images:
    if(exist('h0'))
        if(sum(ismember(findall(0,'type','figure'),h0)))
            close(h0)
        end
    end
    if(exist('h1'))
        if(sum(ismember(findall(0,'type','figure'),h1)))
            close(h1)
        end
    end
    
    disp(['Current Roi Image, Num Rois: ' num2str(roi_data.num_rois)]); 
    h0 = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);
    imagesc(roi_data.im_roi); 
    colormap bone
    caxis([-0 nanmedian(nanmedian(roi_data.im_roi(:)))*20])
    axis square
    title(['Num ROIs added: ' num2str(roi_data.num_rois) '  Add ROI? y/n']); 
    in = input('Want to Add ROI? y/n:   ', 's');
    in = lower(in);
    if(isempty(in) || strcmp(in, 'y'))
        disp('Draw the ROI (click, hold click, draw on image)...')
        title('Draw ROI'); 
        numArea = 1; %Only draw one ROI
        
        draw_complete = 0;
        while(~draw_complete)
            drawRois(numArea);  %draw rois

            %Collect roi data: 
            rois = [];
            hF = gcf;
            hP = findobj(hF, 'Tag', 'ROIPatch'); %finds rois
            for rr = 1:size(hP, 1)
                rois(:,:,rr) = getUD(hP(rr), 'binroi'); %gets rois
            end
            if(~isempty(rois))
                draw_complete = 1; 
            end
        end
        
        %Add single ROI 
        Im_roi_i = roi_data.im_bg;
        Im_roi_i(:,:,3) = rois; 
        h1 = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);
        imagesc(Im_roi_i); 
        colormap bone
        caxis([-0 nanmedian(nanmedian(roi_data.im_roi(:)))*20])
        axis square; title(['Added Candidate Roi # ' num2str(roi_data.num_rois+1) ' Keep? y/n']);
        %Check whether to keep ROI:
        in = lower(input('Keep ROI? y/n:   ', 's'));
        if(isempty(in) || strcmp(in, 'y'))
            %Update general information: 
            roi_data.num_rois = roi_data.num_rois+1; 
            %This ROI mask binary:
            roi_data.roi_bin_cell{roi_data.num_rois}      = rois;
            %Update x,y,r: 
            [roi_ctr] = roi_bin_cell2center_radius({rois});
            roi_data.x = [roi_data.x roi_ctr.x];
            roi_data.y = [roi_data.y roi_ctr.y];
            roi_data.r = [roi_data.r roi_ctr.r]; 
            
            roi_idxs = find(rois);
            roi_data.roi_mask(roi_idxs)        = roi_data.num_rois;
            roi_data.roi_mask_bin(roi_idxs)    = 1;   
            roi_data.im_roi(:,:,3)               = roi_data.roi_mask_bin;
            
            %Update channel information: 
            chan_selected_bool = 0;
            while(~chan_selected_bool)
                chan_vec = [0; 1]; 
                roi_data.chan_logical = [roi_data.chan_logical chan_vec];                     
                g_mod = squeeze(roi_data.im_roi_rg(:,:,2));
                g_mod(roi_idxs) = 1;
                roi_data.im_roi_rg(:,:,2) = g_mod;
                chan_selected_bool = 1; 
            end            
        end
    else
        roi_complete_bool = 1;
        disp('done')
        if(sum(ismember(findall(0,'type','figure'),h0)))
            close(h0)
        end
        if(sum(ismember(findall(0,'type','figure'),h1)))
            close(h1)
        end
                
        %Channel specific information: 
        roi_data = roi_data2chan(roi_data); 
        
        h = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);
        imagesc(roi_data.im_roi_rg); 
        colormap bone
        caxis([-0 nanmedian(nanmedian(roi_data.im_roi(:)))*20])
        axis square
        title(['ROI addition complete! Num ROIs: ' num2str(roi_data.num_rois)]);  
        
        h = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);
        imagesc(roi_data.im_roi);
        colormap bone
        caxis([-0 1000])
        axis square
        title(['ROI addition complete! Num ROIs: ' num2str(roi_data.num_rois)]);         
    end
end

end