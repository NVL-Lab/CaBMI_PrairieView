function bmi_raw_data = clean_bmi_raw_data(bmi_raw_data_file)
%{
Function to remove the last nans values of bmi_raw_data
%}

    % load data
    bmi_raw_data = load(bmi_raw_data_file);
    % Get the fieldnames of bmi_raw_data.data
    data_fields = fieldnames(bmi_raw_data.data);

    % Iterate through the fields
    for i = 1:numel(data_fields)
        % Check if the field contains a 2D array
        if ismatrix(bmi_raw_data.data.(data_fields{i}))
            % Get the size of the array
            array_size = size(bmi_raw_data.data.(data_fields{i}));

            % Check if either dimension N or M is larger than bmi_raw_data.frame
            if array_size(1) > bmi_raw_data.data.frame || array_size(2) > bmi_raw_data.data.frame
                % Determine the number of components in the M dimension to remove
                components_to_remove = array_size(2) - bmi_raw_data.data.frame;

                % Remove the components from the M dimension
                bmi_raw_data.data.(data_fields{i}) = ...
                    bmi_raw_data.data.(data_fields{i})(:, 1:end-components_to_remove);
            end
        end
    end