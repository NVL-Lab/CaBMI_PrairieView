function [roi_data_sel, sel_idxs] = select_roi_data(roi_data, sel_idxs)
%Example roi_data:
% roi_data = 
% 
%   struct with fields:
% 
%         num_chan: 1
%            im_bg: [512×512×3 double]
%         num_rows: 512
%         num_cols: 512
%         num_rois: 43
%         roi_mask: [512×512 double]
%     roi_mask_bin: [512×512 double]
%     roi_bin_cell: {1×43 cell}
%     chan_logical: [2×43 double]
%                x: [1×43 double]
%                y: [1×43 double]
%                r: [1×43 double]
%           im_roi: [512×512×3 double]
%        im_roi_rg: [512×512×3 double]
%             chan: [1×2 struct]

% sel_idxs                        = unique(sel_idxs);
roi_data_sel.sel_idxs           = sel_idxs; %ascending order
roi_data_sel.im_bg              = roi_data.im_bg; 
roi_data_sel.num_rows           = roi_data.num_rows; 
roi_data_sel.num_cols           = roi_data.num_cols; 
roi_data_sel.num_rois           = length(sel_idxs); 
roi_data_sel.roi_bin_cell       = roi_data.roi_bin_cell(sel_idxs); 
roi_data_sel.chan_logical       = roi_data.chan_logical(:,sel_idxs);
roi_data_sel.x                  = roi_data.x(sel_idxs);
roi_data_sel.y                  = roi_data.y(sel_idxs);
roi_data_sel.r                  = roi_data.r(sel_idxs);


roi_data_sel.im_roi     = roi_data.im_bg; 
roi_data_sel.im_roi_rg  = roi_data.im_bg; 
roi_data_sel.roi_mask   = zeros(roi_data_sel.num_rows, roi_data_sel.num_cols); 
roi_data_sel.roi_mask_bin   = zeros(roi_data_sel.num_rows, roi_data_sel.num_cols); 
    
for i = 1:roi_data_sel.num_rois
    roi_num = sel_idxs(i);
    roi_i = roi_data_sel.roi_bin_cell{i};
    roi_idxs = find(roi_i);
    roi_data_sel.roi_mask(roi_idxs)        = roi_num;
    roi_data_sel.roi_mask_bin(roi_idxs)    = 1;   
        
	%Update im_roi_rg: 
	%         del_data.chan_logical
	chan_idx = find(roi_data_sel.chan_logical(:, i)); 
	%         chan_idx
	chan_im = squeeze(roi_data_sel.im_roi_rg(:,:,chan_idx)); 
	chan_im(roi_idxs) = 1; 
	roi_data_sel.im_roi_rg(:,:,chan_idx) = chan_im;  
end
roi_data_sel.im_roi(:,:,3)               = roi_data_sel.roi_mask_bin; 

%Visualize:
screen_size = get(0,'ScreenSize');
h = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);
imagesc(roi_data_sel.im_roi); axis square
title('Background Image + Rois'); 
%
h = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);
imagesc(roi_data_sel.roi_mask); axis square
title(['Num ROI: ' num2str(roi_data_sel.num_rois)]); 
    
 