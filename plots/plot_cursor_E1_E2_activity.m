function [h, offset_vec] = plot_cursor_E1_E2_activity(cursor, E1, E2, n, E_id, E_color, offset)
%4.18.19
%n - neural activity, num_samples x num_neurons
%E_id - vector of 1's and 2's, indicating membership to E1 or E2
%E_color - cell array of colors for E1 and E2

%E1, E2 activity:
num_E1 = sum(E_id==1); 
E1_sel = find(E_id == 1); 
n_E1 = n(:,E1_sel); 

num_E2 = sum(E_id==2); 
E2_sel = find(E_id == 2); 
n_E2 = n(:,E2_sel); 

%-------------------------------
if nargin < 7
    offset = 0;
end

offset_vec = [offset]; 
%CURSOR
h = figure; hold on;
plot(cursor); 
y_amp = max(cursor)-min(cursor); 
offset = offset + y_amp; 
offset_vec = [offset_vec offset]; 

%E1
E1_plot = E1-offset;
plot(E1_plot, 'k'); 
y_amp = max(E1_plot)-min(E1_plot); 
offset = offset + y_amp; 
offset_vec = [offset_vec offset]; 

%E2
E2_plot = E2-offset;
plot(E2_plot, 'r'); 
y_amp = max(E2_plot)-min(E2_plot); 
offset = offset + y_amp; 
offset_vec = [offset_vec offset]; 
legend({'cursor', 'E1', 'E2'}); 

%Plot E1
for i=1:num_E1
    y_plot = n_E1(:,i); 
    y_plot = y_plot-min(y_plot); 
    y_amp = max(y_plot); 
    
    offset = offset + y_amp;
    offset_vec = [offset_vec offset]; 
    y_plot = y_plot-offset;
    plot(y_plot, 'Color', E_color{1});     
end

%Plot E2
for i=1:num_E2
    y_plot = n_E2(:,i); 
    y_plot = y_plot-min(y_plot); 
    y_amp = max(y_plot); 
    
    offset = offset + y_amp;
    offset_vec = [offset_vec offset]; 
    y_plot = y_plot-offset;
    plot(y_plot, 'Color', E_color{2});     
end

% num_neurons = size(n,2); 
% if nargin < 4
%     offset = 0;
% end
% offset_vec = [offset]; 
% h = figure; hold on;
% for i=1:num_neurons
%     y_plot = n(:,i);
%     y_plot = y_plot-min(y_plot); 
%     y_amp = max(y_plot); 
%     if(i>1)
%         offset = offset + y_amp;
%         offset_vec = [offset_vec offset]; 
%     end
%     y_plot = y_plot-offset;
%     plot(y_plot, 'Color', E_color{E_id(i)}); 
% end