function data_path = find_data_path(folder_list, component)
%{ 
Function that creates the correct folder_path after obtaining which 
 location in the computer each hard-drive has (given by folder_list).
%}

    %% folder_paths structure
    % hardcoded structure given the location of the data in the same HD
    folder_paths = struct();
    folder_paths.FA = {'m13', 'm15', 'm16', 'm18', 'm25'};
    folder_paths.FB = {'m21', 'm22', 'm26'};
    folder_paths.FC = {'m23', 'm27', 'm28', 'm29'};

    % Find the corresponding key in folder_list
    folder_key_list = fieldnames(folder_list);
    matching_key = '';
    for i = 1:numel(folder_key_list)
        data_cell = folder_paths.(folder_key_list{i});
        if any(strcmp(data_cell, component))
            matching_key = folder_key_list{i};
            break;
        end
    end

    % Check if a matching key was found
    if isempty(matching_key)
        data_path = '';
        disp('No matching key found in folder_list.');
        return;
    end

   
    % Create the full data path using folder_list value and matching component
    data_path = folder_list.(matching_key);
end