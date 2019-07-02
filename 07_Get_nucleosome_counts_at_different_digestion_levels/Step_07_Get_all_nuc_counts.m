sample_label = {'1min', '2min', '5min', '15min', '40min', '60min'};
noSamples = numel(sample_label);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the normalized counts %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% S2_exp1
pathName = '../../data/S2_exp1';
load(fullfile(pathName, 'Typical_nucleosome_positions_in_S2_exp1.mat'), 'Chr', 'Loc', ...
    'Occ_score', 'Linker_score', 'position_std')
load(fullfile(pathName, 'local_norm_factor.S2_exp1.mat'), 'local_norm_factor')
aligned_local_norm_factor = align_sites(local_norm_factor, Chr, Loc, 50, 50);
local_norm_factor_per_nuc = median(aligned_local_norm_factor, 2);

% Remove nucs. from bad regions
badNucs = isnan(local_norm_factor_per_nuc);

Chr(badNucs) = [];
Loc(badNucs) = [];
Occ_score(badNucs) = [];
Linker_score(badNucs) = [];
position_std(badNucs) = [];
noNucs = numel(Loc);

aligned_local_norm_factor = align_sites(local_norm_factor, Chr, Loc, 50, 50);
local_norm_factor_per_nuc = median(aligned_local_norm_factor, 2);

NucCount = nan(noNucs, noSamples);
for s = 1:noSamples
    load(fullfile(pathName, ['Spikein_normalized_profiles.100-200.S2_exp1_', sample_label{s}, '.dm6.mat']), 'Occ_147')
    
    alignedOcc = align_sites(Occ_147, Chr, Loc, 50, 50);
    NucCount(:, s) = median(alignedOcc, 2) ./ local_norm_factor_per_nuc;
end

% Save counts
save(fullfile(pathName, 'NucCounts.100-200.S2_exp1.mat'), 'Chr', 'Loc', 'NucCount', 'noNucs', ...
    'Occ_score', 'Linker_score', 'position_std');



%% S2_exp2
pathName = '../../data/S2_exp2';
load(fullfile(pathName, 'Typical_nucleosome_positions_in_S2_exp2.mat'), 'Chr', 'Loc', ...
    'Occ_score', 'Linker_score', 'position_std')
load(fullfile(pathName, 'local_norm_factor.S2_exp2.mat'), 'local_norm_factor')
aligned_local_norm_factor = align_sites(local_norm_factor, Chr, Loc, 50, 50);
local_norm_factor_per_nuc = median(aligned_local_norm_factor, 2);

% Remove nucs. from bad regions
badNucs = isnan(local_norm_factor_per_nuc);

Chr(badNucs) = [];
Loc(badNucs) = [];
Occ_score(badNucs) = [];
Linker_score(badNucs) = [];
position_std(badNucs) = [];
noNucs = numel(Loc);

aligned_local_norm_factor = align_sites(local_norm_factor, Chr, Loc, 50, 50);
local_norm_factor_per_nuc = median(aligned_local_norm_factor, 2);

NucCount = nan(noNucs, noSamples);
for s = 1:noSamples
    load(fullfile(pathName, ['Spikein_normalized_profiles.100-200.S2_exp2_seq1_', sample_label{s}, '.dm6.mat']), 'Occ_147')
    
    alignedOcc = align_sites(Occ_147, Chr, Loc, 50, 50);
    NucCount(:, s) = median(alignedOcc, 2) ./ local_norm_factor_per_nuc;
end

% Save counts
save(fullfile(pathName, 'NucCounts.100-200.S2_exp2.mat'), 'Chr', 'Loc', 'NucCount', 'noNucs', ...
    'Occ_score', 'Linker_score', 'position_std');



%% Kc167_exp1
pathName = '../../data/Kc167_exp1';
load(fullfile(pathName, 'Typical_nucleosome_positions_in_Kc167_exp1.mat'), 'Chr', 'Loc', ...
    'Occ_score', 'Linker_score', 'position_std')
load(fullfile(pathName, 'local_norm_factor.Kc167_exp1.mat'), 'local_norm_factor')
aligned_local_norm_factor = align_sites(local_norm_factor, Chr, Loc, 50, 50);
local_norm_factor_per_nuc = median(aligned_local_norm_factor, 2);

% Remove nucs. from bad regions
badNucs = isnan(local_norm_factor_per_nuc);

Chr(badNucs) = [];
Loc(badNucs) = [];
Occ_score(badNucs) = [];
Linker_score(badNucs) = [];
position_std(badNucs) = [];
noNucs = numel(Loc);

aligned_local_norm_factor = align_sites(local_norm_factor, Chr, Loc, 50, 50);
local_norm_factor_per_nuc = median(aligned_local_norm_factor, 2);

NucCount = nan(noNucs, noSamples);
for s = 1:noSamples
    load(fullfile(pathName, ['Spikein_normalized_profiles.100-200.Kc167_exp1_seq1_', sample_label{s}, '.dm6.mat']), 'Occ_147')
    
    alignedOcc = align_sites(Occ_147, Chr, Loc, 50, 50);
    NucCount(:, s) = median(alignedOcc, 2) ./ local_norm_factor_per_nuc;
end

% Save counts
save(fullfile(pathName, 'NucCounts.100-200.Kc167_exp1.mat'), 'Chr', 'Loc', 'NucCount', 'noNucs', ...
    'Occ_score', 'Linker_score', 'position_std');


