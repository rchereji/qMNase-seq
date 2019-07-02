function RescaledProfile = local_normalization_cell_array_to_avg_1(Profile)

chrLen = cellfun(@length, Profile);

% Initialize an empty cell array
RescaledProfile = cell(size(Profile));

% Rescale the cell array by chromosome averages
noChr = numel(Profile);
for chr = 1:noChr
    % Compute the chromosome average (discarding regions with artifacts - very high coverage)
    tmp = Profile{chr};
    threshold = 5; % 5 times the genome average
    
    aboveThreshold = (tmp > threshold);  %where above threshold
    % aboveThreshold is a logical array, where 1 when above threshold, 0, below.
    % we thus want to calculate the difference between rising and falling edges
    aboveThreshold = [false, aboveThreshold, false];  %pad with 0's at ends
    edges = diff(aboveThreshold);
    rising = find(edges == 1);     % rising/falling edges
    falling = find(edges == -1);
    spanWidth = falling - rising;  % width of span of 1's (above threshold)
    wideEnough = spanWidth >= 5;
    startPos = max(rising(wideEnough) - 50, 1);             % start of each span - 50bp
    endPos = min(falling(wideEnough)-1 + 50, chrLen(chr));  % end of each span + 50bp
    
    % Normalize each regions individually
    noRegions = numel(startPos);
    for r = 1:noRegions
        tmp(startPos(r):endPos(r)) = tmp(startPos(r):endPos(r)) / nanmean(tmp(startPos(r):endPos(r)));
    end
    
    mov_average_2k = movmean(tmp, 2001);
    
    RescaledProfile{chr} = tmp ./ mov_average_2k;
end
