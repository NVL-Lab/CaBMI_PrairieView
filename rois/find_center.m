function [x,y]=find_center(temp, mat, toplot)
%{
Function to find the center of mass for the detected cells
temp is the template / roi_mask 
mat is the raw image
toplot is a flag to allow plotting
%}
    if nargin < 2
        toplot = false;
        mat = 0;
    elseif nargin < 3
        toplot = true;
    end
    cell_ind = unique(temp(:));
    cell_ind(cell_ind==0) = []; 
    num_cells = length(cell_ind);
    y=zeros(num_cells,1);
    x=zeros(num_cells,1);
    for i=1:num_cells
        y(i)= round(mean(find(temp'==cell_ind(i)))/size(temp,2));
        x(i)= round(mean(find(temp ==cell_ind(i)))/size(temp,1));
    end

    if toplot
        figure ,
        if ~isempty(y)
            subplot (1,2,1)
            imagesc(temp), colormap bone, caxis([-0 nanmedian(nanmedian(temp(:)))*5]), hold on, scatter (x,y, 'filled', 'r'), hold off, axis square            
%             imagesc(temp), colormap bone, caxis([0 1]), hold on, scatter (x,fliplr(y), 'filled', 'r'), hold off, axis square

            subplot (1,2,2)
            imagesc(mat), colormap bone, caxis([-0 nanmedian(nanmedian(mat(:)))*5]), hold on, scatter (x,y, 'filled', 'r'), hold off, axis square            
%             imagesc(mat), colormap bone, caxis([-0 nanmean(nanmean(mat))*4]), hold on, scatter (x,fliplr(y), 'filled', 'r'), hold off, axis square
        else
            subplot (1,2,1)
            imagesc(temp), colormap bone, caxis([-0 nanmedian(nanmedian(temp(:)))*5])

            subplot (1,2,2)
            imagesc(mat), colormap bone, caxis([-0 nanmedian(nanmedian(mat(:)))*5])

        end
    end
end