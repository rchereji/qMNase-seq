function RescaledProfile = normalize_cell_array_to_avg_1(Profile)

% Initialize an empty cell array
RescaledProfile = cell(size(Profile));

% Rescale the cell array by chromosome averages
noChr = numel(Profile);
for chr = 1:noChr
    % Compute the chromosome average (discarding regions with artifacts - very high coverage)
    tmp = Profile{chr};
    Filter = tmp > 100 * nanmean(tmp);
    Filter = imdilate(Filter, strel('line', 101, 0));
    tmp(Filter) = nan;
    chrAverage = nanmean(tmp);
    
    % Rescale profile by chromosome average
    RescaledProfile{chr} = Profile{chr} / chrAverage;
end