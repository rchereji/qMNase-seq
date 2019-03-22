sample_label = {'1min', '2min', '5min', '15min', '40min', '60min'};
noSamples = numel(sample_label);

load('data/Typical_nucleosome_positions_in_dm6.mat', 'Chr', 'Loc')
load('data/local_norm_factor.dm6.mat', 'local_norm_factor')
aligned_local_norm_factor = align_dm6_sites(local_norm_factor, Chr, Loc, 50, 50);
local_norm_factor_per_nuc = median(aligned_local_norm_factor, 2);

NucCount = nan(length(Loc), noSamples);
for s = 1:noSamples
    load(['data/Spikein_normalized_profiles.100-200.S2_', sample_label{s}, '.dm6.mat'], 'Occ_147')
    
    alignedOcc = align_dm6_sites(Occ_147, Chr, Loc, 50, 50);
    NucCount(:, s) = median(alignedOcc, 2) ./ local_norm_factor_per_nuc;
end

%% Remove nucs. from bad regions
badNucs = isnan(local_norm_factor_per_nuc);

NucCount(badNucs, :) = [];
Chr(badNucs) = [];
Loc(badNucs) = [];
noNucs = numel(Loc);

save('data/NucCounts.mat', 'Chr', 'Loc', 'NucCount', 'noNucs');
