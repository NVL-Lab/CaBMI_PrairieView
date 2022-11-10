function [roi_data] = delete_roi_2chan(plot_images, roi_data)
%Enter idxs you want deleted: 

screen_size = get(0,'ScreenSize');
%Show im_roi_rg, and roi_mask

% Show the red and green channels:
close all;
for plot_i=1:length(plot_images)
    im_plot = plot_images(plot_i).im; 
    im_title = plot_images(plot_i).label; 
    h = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);
    imagesc(im_plot); axis square; colormap('gray'); title(im_title);    
end

%Ask for user input on idxs to delete
%show updated roi
%ask if they want to accept
%ask if they are done deleting

complete_bool = 0; 
while(~complete_bool)
    %Show images
    h1 = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);
    imagesc(roi_data.im_roi_rg); axis square; title('roi colored by rg'); 
    h2 = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);
    imagesc(roi_data.roi_mask); axis square; title('roi idxs'); 
    
    del_idxs = input('enter roi idxs to delete in vector [roi_1 roi_2 ... ]: ');
    del_data = delete_rois(roi_data, del_idxs);

    %Show new images: 
    h3 = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);
    imagesc(del_data.im_roi_rg); axis square; title('roi colored by rg, AFTER DELETION'); 
    h4 = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);
    imagesc(del_data.roi_mask); axis square; title('roi idxs, AFTER DELETION')    
    
    in = input('UNDO? y/n:   ', 's');
    in = lower(in);
    if(isempty(in) || strcmp(in, 'n'))      
        roi_data = del_data; 
    end
    
    in = input('Done Deleting? y/n:   ', 's');
    in = lower(in);
    if(isempty(in) || strcmp(in, 'y'))      
        complete_bool = 1; 
    end
    
    %Close images: 
    %Close images:
    if(exist('h1'))
        if(sum(ismember(findall(0,'type','figure'),h1)))
            close(h1)
        end
    end
    if(exist('h2'))
        if(sum(ismember(findall(0,'type','figure'),h2)))
            close(h2)
        end
    end
    %Close images:
    if(exist('h3'))
        if(sum(ismember(findall(0,'type','figure'),h3)))
            close(h3)
        end
    end
    if(exist('h4'))
        if(sum(ismember(findall(0,'type','figure'),h4)))
            close(h4)
        end
    end    
end

end


function del_data = delete_rois(roi_data, del_idxs)
        %1) remove from roi_bin_cell
        %2) remove from idx lists, decrement roi count
        %3) rebuild im_roi, im_roi_rg, roi_mask, roi_mask_bin: 
        %3a) pool
        %3b) r
        %3c) g
    
    del_data = roi_data;
    del_data.roi_bin_cell(del_idxs)         = []; 
    del_data.chan_logical(:,del_idxs)       = []; 
    del_data.x(del_idxs)                    = []; 
    del_data.y(del_idxs)                    = []; 
    del_data.r(del_idxs)                    = []; 
    del_data.num_rois                       = length(del_data.roi_bin_cell); 
    
    %Rebuild: 
    %im_roi, im_roi_rg, roi_mask, roi_mask_bin
    del_data.im_roi     = del_data.im_bg; 
    del_data.im_roi_rg  = del_data.im_bg; 
    del_data.roi_mask   = zeros(del_data.num_rows, del_data.num_cols); 
    del_data.roi_mask_bin   = zeros(del_data.num_rows, del_data.num_cols); 
    
    for i = 1:del_data.num_rois
        roi_i = del_data.roi_bin_cell{i};
%         size(roi_i)
        roi_idxs = find(roi_i);
        
        del_data.roi_mask(roi_idxs)        = i;
        del_data.roi_mask_bin(roi_idxs)    = 1;   
        
        %Update im_roi_rg: 
%         del_data.chan_logical
        chan_idx = find(del_data.chan_logical(:, i)); 
%         chan_idx
        chan_im = squeeze(del_data.im_roi_rg(:,:,chan_idx)); 
        chan_im(roi_idxs) = 1; 
        del_data.im_roi_rg(:,:,chan_idx) = chan_im;  
    end
    del_data.im_roi(:,:,3)               = del_data.roi_mask_bin; 
    
    %Update channel information: 
    [del_data] = roi_data2chan(del_data);
    
end

%                     roi_data.chan_logical = [roi_data.chan_logical chan_vec]; 
%                     r_mod = squeeze(roi_data.im_roi_rg(:,:,1));
%                     r_mod(roi_idxs) = 1;
%                     roi_data.im_roi_rg(:,:,1) = r_mod;                  
%                     chan_selected_bool = 1;