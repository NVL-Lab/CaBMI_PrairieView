function roi_data = label_mask2roi_data_single_channel(im_bg, label_mask, chan_data)
%Currently supports only num_chan = 1
%im_bg: num_row X num_col X 1
%TODO: generalize to multi channel
%INPUT:
%label_mask - comes from automatic ROI detection. 
%EG:
%    [mask_intermediate, ~] = imFindCellsTM (im_bg, template_diam, thres, cell_diam, finemode, temmode);
%     init_roi_mask = bwlabel(mask_intermediate);
% matrix of size(im_bg), with ROI's going from 1 to num_roi's

screen_size = get(0,'ScreenSize');
num_chan = length(chan_data);

roi_data.num_chan       = num_chan;
roi_data.im_bg          = repmat(im_bg, [1 1 3]); 
roi_data.num_rows       = size(im_bg,1);
roi_data.num_cols       = size(im_bg,2);

roi_data.num_rois       = max(label_mask(:)); 
roi_data.roi_mask       = label_mask; 
roi_data.roi_mask_bin   = label_mask > 0; 
roi_data.roi_bin_cell   = {};
roi_data.chan_logical   = ...
    [zeros(1,roi_data.num_rois); ones(1,roi_data.num_rois)]; %num_chan x num_roi

for roi_i = 1:roi_data.num_rois
    roi_data.roi_bin_cell{roi_i} = ...
        label_mask == roi_i; 
end

%Find x,y,r: 
[roi_ctr] = roi_bin_cell2center_radius(roi_data.roi_bin_cell);
roi_data.x              = roi_ctr.x; 
roi_data.y              = roi_ctr.y; 
roi_data.r              = roi_ctr.r; 

%im_roi
roi_data.im_roi         = zeros(roi_data.num_rows, roi_data.num_cols, 3); 
roi_data.im_roi(:,:,2)  = im_bg; 
roi_data.im_roi(:,:,3)  = roi_data.roi_mask_bin;

%im_roi_rg
roi_idxs                = find(roi_data.roi_mask_bin); 
roi_data.im_roi_rg      = repmat(im_bg, [1 1 3]); 
g_mod                   = squeeze(roi_data.im_roi_rg(:,:,2));
g_mod(roi_idxs)         = 1;
roi_data.im_roi_rg(:,:,2) = g_mod;


for i = 1:num_chan
    chan_i                                  = chan_data(i).chan_idx; 
    roi_data.chan(chan_i).label             = chan_data(i).label;
    roi_data.chan(chan_i).num_rois          = roi_data.num_rois; 
    roi_data.chan(chan_i).idxs              = 1:roi_data.num_rois; 
    roi_data.chan(chan_i).im_roi            = roi_data.im_roi; 
    roi_data.chan(chan_i).roi_mask          = roi_data.roi_mask;     
    roi_data.chan(chan_i).roi_mask_bin      = roi_data.roi_mask_bin;     
end

%
h = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);
imagesc(roi_data.im_bg), caxis([-0 nanmean(nanmean(im_bg(:)))*2]); axis square
title('Background Image'); 
%
h = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);
imagesc(roi_data.im_roi), caxis([-0 nanmean(nanmean(im_bg(:)))*2]); axis square
title(['Num ROI: ' num2str(roi_data.num_rois)]); 


