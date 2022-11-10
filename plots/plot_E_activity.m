function [h, offset_vec] = plot_E_activity(t,n, E_id, E_color, offset)
%4.18.19
%n - neural activity, num_samples x num_neurons
%E_id - vector of 1's and 2's, indicating membership to E1 or E2
%E_color - cell array of colors for E1 and E2
num_neurons = size(n,2); 
if nargin < 5
    offset = 0;
end
offset_vec = [offset]; 
h = figure; hold on;
for i=1:num_neurons
    y_plot = n(:,i);
    y_plot = y_plot-min(y_plot); 
    y_amp = max(y_plot); 
    if(i>1)
        offset = offset + y_amp;
        offset_vec = [offset_vec offset]; 
    end
    y_plot = y_plot-offset;
%     size(t)
%     size(y_plot)
    plot(t,y_plot, 'Color', E_color{E_id(i)}); 
end