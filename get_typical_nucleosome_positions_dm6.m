% This script searches for the typical positions of well positioned
% nucleosomes in Drosophila, by identifying the medians of all nucleosome
% dyad clusters

%% Load dyad data
load('data/All_Dyads.100-200.S2_cells.dm6.mat', 'Dyads', 'chrLen')
Occ_101 = extend_dyads_to_fixed_footprint(Dyads, 101);
Occ_101 = normalize_cell_array_to_avg_1(Occ_101);

%% Initialize the positions to be equally spaced (50 bp apart) and to cover the whole chromosome 
Chr = [];
Loc = [];

noChr = numel(chrLen);
for c = 1:noChr
    tmpLoc = [50:50:chrLen(c)]';
    Chr = [Chr; c*ones(numel(tmpLoc), 1)];
    Loc = [Loc; tmpLoc];
end

%% Search for the medians of dyad clusters; iterate the search process until the shift in medians is at most 1 bp
noSites = numel(Chr);
nucShift = 100 * ones(noSites, 1); % Initialize with a large number
maxShift = round(max(abs(nucShift)));

noIter = 0;
while maxShift > 1
    AlignedDyads = align_dm6_sites(Dyads, Chr, Loc, 73, 73);
    AlignedDyads(isnan(AlignedDyads)) = 0;
    
    % Find medians and update the nuc positions to the new medians
    x = -73:73;
    for s = 1:noSites
        u = repelem(x, AlignedDyads(s,:));
        nucShift(s) = round(median(u));
        Loc(s) = Loc(s) + nucShift(s);
    end
    
    % Get the maximum nuc. shift
    maxShift = round(max(abs(nucShift)));
    noIter = noIter + 1;
    if maxShift > 1
        fprintf('Iteration %d done. Max. shift: %d. Continue.\n', noIter, maxShift)
    else
        fprintf('Iteration %d done. Stop.\n', noIter)
    end
end

%% Keep only the unique positions
M = unique([Chr, Loc], 'rows');
distBwNucs = [0; diff(M(:,2))];
rowsToEliminate = (isnan(M(:,2))) | (M(:,2) < 73) | (distBwNucs < 5);
M = M(~rowsToEliminate, :);

Chr = M(:,1);
Loc = M(:,2);

%% Compute average occupancy and linker scores
AlignedOcc = align_dm6_sites(Occ_101, Chr, Loc, 73+10, 73+10);

noSites = numel(Chr);
OccScore = nan(noSites, 1);
LinkerScore = nan(noSites, 1);
for s = 1:noSites
    OccScore(s) = nanmean(AlignedOcc(s, 73+10+[-50:50]));
    LinkerScore(s) = nanmean(AlignedOcc(s, [1:10, end-9:end])) ./ nanmean(AlignedOcc(s, 73+10+[-50:50]));
end

%%
% Eliminate peak calls with low occupancy scores or high linker scores
badCalls = (OccScore < 0.15) | (LinkerScore > 1);
Chr = Chr(~badCalls);
Loc = Loc(~badCalls);
OccScore = OccScore(~badCalls);
LinkerScore = LinkerScore(~badCalls);

%% Get the clustered peak calls and keep only the one with the highest occupancy score
distBwNucs = [NaN; diff(Loc)];
distBwNucs(distBwNucs < 0) = NaN;

nucsInClusters = (distBwNucs < 75);
lastNucInClusters = strfind(nucsInClusters', [1,0]) + 1;
nucsInClusters(lastNucInClusters) = 1;
[Label, noClust] = bwlabel(nucsInClusters);

noSites = numel(Chr);
badCalls = zeros(noSites, 1);

for c = 1:noClust
    NucInd = (Label == c);
    Scores = OccScore(NucInd);
    [~, Ind] = max(Scores);
    NucInd(Ind) = [];
    badCalls(NucInd) = 1;
end
    
Chr = Chr(~badCalls);
Loc = Loc(~badCalls);
OccScore = OccScore(~badCalls);
LinkerScore = LinkerScore(~badCalls);

IGV_Score = round(1000 * (1-LinkerScore));

%% Save the typical nucleosome positions in a MATLAB file
save('data/Typical_nucleosome_positions_in_dm6.mat', 'Chr', 'Loc', 'OccScore', 'LinkerScore', 'IGV_Score')

%% Save the typical nucleosome positions in a BED file
load('data/chrLen_dm6.mat', 'chrReference')
ChrLabel = chrReference(Chr);

noSites = numel(Chr);
BEDcontents = 'track name=typicalNucPos description="Typical nucleosome positions" useScore=1 itemRgb="Off"\n';
for s = 1:noSites
    BEDcontents = [BEDcontents, sprintf('%s\t%d\t%d\tNuc_%d\t%d\t*\t%d\t%d\t0,0,255\n', ...
        ChrLabel{s}, Loc(s)-73, Loc(s)+73, s, IGV_Score(s), Loc(s)-73, Loc(s)+73)];
end

% Save BED file
fileID = fopen('data/Typical_nucleosome_positions_in_dm6.bed', 'w');
fprintf(fileID, '%s', BEDcontents);
fclose(fileID);
