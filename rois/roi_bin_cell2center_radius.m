function [roi_ctr] = roi_bin_cell2center_radius(roi_bin_cell)
num_roi = length(roi_bin_cell); 
roi_ctr.x = zeros(1, num_roi); 
roi_ctr.y = zeros(1, num_roi); 
roi_ctr.r = zeros(1, num_roi); 

for roi_i = 1:num_roi
    roi_im = roi_bin_cell{roi_i};
    roi_ctr.x(roi_i)    = round(mean(find(roi_im))/size(roi_im,1));
    roi_ctr.y(roi_i)    = round(mean(find(roi_im'))/size(roi_im,2));
    x_occupy = find(sum(roi_im, 1)>0); 
    y_occupy = find(sum(roi_im, 2)>0); 
    xdel = max(x_occupy)-min(x_occupy);
    ydel = max(y_occupy) - min(y_occupy);
    roi_ctr.r(roi_i)    = max(xdel, ydel)/2; 
end