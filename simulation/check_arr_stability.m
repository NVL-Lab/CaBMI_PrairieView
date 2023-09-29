function isStable = check_arr_stability(arr, num_samp, threshold)
    %{
    Function to check the stability of a neuron
    %}    

    if nargin <3
        threshold = 0.2;
    end
    if nargin < 2
        num_samp = 10000;
    end

    if length(arr) < num_samp
        isStable = true;  % Not enough data points to calculate stability.
        return;
    end

    indices = 1 : num_samp/2 : length(arr);
    std_arr = zeros(1, length(indices) - 2);

    for i = 2 : length(indices) - 1
        startIndex = indices(i) - num_samp/2;
        endIndex = indices(i) + num_samp/2;
        std_arr(i - 1) = nanstd(arr(startIndex:endIndex));
    end

    mean_max = max([nanmean(std_arr) * (1 + threshold), nanmean(std_arr) + 0.2]);
    mean_min = min([nanmean(std_arr) * (1 - threshold), nanmean(std_arr) - 0.2]);

    stability = sum(std_arr > mean_max) + sum(std_arr < mean_min);

    if stability > 0
        isStable = false;
    else
        isStable = true;
    end
end
