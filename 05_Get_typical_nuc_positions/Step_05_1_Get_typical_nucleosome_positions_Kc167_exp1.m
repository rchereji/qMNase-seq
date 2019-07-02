% This script searches for the typical positions of well positioned
% nucleosomes in Drosophila, by identifying the medians of all nucleosome
% dyad clusters

% Multiple options are possible: median, mean, mode
% We will search for the median positions for each cluster of nuc. centers
cluster_center_type = 'median';
cluster_center = str2func(cluster_center_type);

expName  = 'Kc167_exp1';
pathName = '../../data/Kc167_exp1';

%% Load dyad data
load(fullfile(pathName, ['All_Raw_Dyads.100-200.', expName, '.dm6.mat']), 'Dyads')
chrLen = cellfun(@length, Dyads);

% Extend dyads to a footprint of 101 bp
Occ_101 = extend_dyads_to_fixed_footprint(Dyads, 101);
% Normalize the coverage for each chromosome
Occ_101 = normalize_cell_array_to_avg_1(Occ_101);
% Do an extra local normalization to correct for local sequencing artifacts
Occ_101 = local_normalization_cell_array_to_avg_1(Occ_101);

%% Initialize the positions to be equally spaced (50 bp apart) and to cover the whole chromosome 
Chr = [];
Loc = [];

noChr = numel(chrLen);
for c = 1:noChr
    tmpLoc = [50:50:chrLen(c)]';
    Chr = [Chr; c*ones(numel(tmpLoc), 1)];
    Loc = [Loc; tmpLoc];
end

%% Search for the median/mean/mode of dyad clusters; 
% iterate the search process until the shift in medians is at most 1 bp

noSites = numel(Chr);
nucShift = 100 * ones(noSites, 1); % Initialize with a large number
maxShift = round(max(abs(nucShift)));

noIter = 0;
while maxShift > 1
    aligned_dyads = align_sites(Dyads, Chr, Loc, 73, 73);
    aligned_dyads(isnan(aligned_dyads)) = 0;
    
    % Find "center" and update the nuc position
    x = -73:73;
    for s = 1:noSites
        u = repelem(x, aligned_dyads(s,:));
        nucShift(s) = round(cluster_center(u));
        Loc(s) = Loc(s) + nucShift(s);
    end
    
    % Get the maximum nuc. shift
    maxShift = round(max(abs(nucShift)));
    noIter = noIter + 1;
    if maxShift > 1
        fprintf('Iteration %d done. Max. shift: %d. Continue...\n', noIter, maxShift)
    else
        fprintf('Iteration %d done. Max. shift: %d. Stop.\n', noIter, maxShift)
    end
end

%% Keep only the unique positions
M = unique([Chr, Loc], 'rows');
distBetweenNucs = [diff(M(:,2)); 0];
rowsToEliminate = (isnan(M(:,2))) | (M(:,2) < 73) | (distBetweenNucs < 5);
M = M(~rowsToEliminate, :);

Chr = M(:,1);
Loc = M(:,2);

%% Compute average occupancy and linker scores
AlignedOcc = align_sites(Occ_101, Chr, Loc, 73+10, 73+10);

noSites = numel(Chr);
Occ_score = nan(noSites, 1);
Linker_score = nan(noSites, 1);
for s = 1:noSites
    Occ_score(s) = nanmean(AlignedOcc(s, 73+10+[-50:50]));
    Linker_score(s) = nanmean(AlignedOcc(s, [1:10, end-9:end])) ./ nanmean(AlignedOcc(s, 73+10+[-50:50]));
end

%% Eliminate peak calls with low occupancy scores or high linker scores
badCalls = (Occ_score < 0.2) | (Linker_score > 1);
Chr = Chr(~badCalls);
Loc = Loc(~badCalls);
Occ_score = Occ_score(~badCalls);
Linker_score = Linker_score(~badCalls);

%% Get the clustered peak calls and keep only the one with the highest occupancy score
distBetweenNucs = [diff(Loc); NaN];
distBetweenNucs(distBetweenNucs < 0) = NaN;

nucsInClusters = (distBetweenNucs < 75);
[Label, noClust] = bwlabel(nucsInClusters);

noSites = numel(Chr);
badCalls = zeros(noSites, 1);

for c = 1:noClust
    NucInd = find(Label == c);
    % Add the last nuc from each cluster
    NucInd = [NucInd; NucInd(end)+1];
    Scores = Occ_score(NucInd);
    [~, Ind] = max(Scores);
    NucInd(Ind) = [];
    badCalls(NucInd) = 1;
end
    
Chr = Chr(~badCalls);
Loc = Loc(~badCalls);
Occ_score = Occ_score(~badCalls);
Linker_score = Linker_score(~badCalls);

IGV_score = round(1000 * (1-Linker_score));

%% Get the standard deviation for each position
aligned_dyads = align_sites(Dyads, Chr, Loc, 73, 73);
aligned_dyads(isnan(aligned_dyads)) = 0;

x = -73:73;
position_std = nan(size(Loc));
noSites = numel(Loc);
for s = 1:noSites
    u = repelem(x, aligned_dyads(s,:));
    position_std(s) = std(u);
end

%% Save results
% Save the typical nucleosome positions in a MATLAB file
save(fullfile(pathName, ['Typical_nucleosome_positions_in_', expName, '.mat']), ...
    'Chr', 'Loc', 'position_std', 'Occ_score', 'Linker_score', 'IGV_score')

% Save the typical nucleosome positions in a BED file
chrReference = {'2L', '2R', '3L', '3R', '4', 'X', 'Y', 'rDNA'};
ChrLabel = chrReference(Chr)';

fileID = fopen(fullfile(pathName, ['Typical_nucleosome_positions_in_', expName, '.bed']), 'W');
header_txt = 'track name="Typical nucleosome positions" description="Typical nucleosome positions" useScore=1 itemRgb="Off"';
fprintf(fileID, '%s\n', header_txt);
Fmt = '%s\t%d\t%d\tNuc_%d\t%d\t*\t%d\t%d\t0,0,255\n';

noSites = numel(Chr);
for s = 1:noSites
  fprintf(fileID, Fmt, ...
      ChrLabel{s}, Loc(s)-74, Loc(s)+73, s, IGV_score(s), Loc(s)-74, Loc(s)+73);
end
fclose(fileID);
