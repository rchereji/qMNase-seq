load('../../data/Kc167_exp1/All_Raw_Dyads.100-200.Kc167_exp1.dm6.mat', 'Dyads')
Occ = extend_dyads_to_fixed_footprint(Dyads, 147);

% Get mean occupancy in sliding 2kb windows 
noChr = numel(Occ);
for chr = 1:noChr
    Occ{chr} = fastsmooth(Occ{chr}, 2001);
end

% Get local normalization factors
genomewide_avg = mean(cat(1, [Occ{:}]));
local_norm_factor = cellfun(@(x) x/genomewide_avg, Occ, 'un', 0);

%%
figure
hist(cat(1, [local_norm_factor{:}]), 0:0.01:10)
xlim([0 10])
% Some regions (2kb windows) have very high/low average occupancy. These
% may represent sequencing artifacts, repeated regions that are not well
% annotated in the genome or missing pieces of the reference genome

%%
% Mark the "bad" regions, and recompute the local normalization factor
% after excluding these regions
bad_region_Filter = cell(size(local_norm_factor));
for chr = 1:noChr
    bad_region_Filter{chr} = (local_norm_factor{chr} < 0.25) | (local_norm_factor{chr} > 4); % 4 times less or more than genomic average
    
    se = strel('line', 301, 0);
    bad_region_Filter{chr} = imdilate(bad_region_Filter{chr}, se);
end

% Check the fraction of the genome that contains problems
mean(cat(1, cat(1, [bad_region_Filter{:}])))
% ans = 0.0372;     
% 3.72% of the genome is marked as bad regions because of extreme coverage 

% Get local normalization factors (after discarding the bad regions)
for chr = 1:noChr
    Occ{chr}(bad_region_Filter{chr}) = NaN;
end

genomewide_avg = nanmean(cat(1, [Occ{:}]));
local_norm_factor = cellfun(@(x) x/genomewide_avg, Occ, 'un', 0);

%%
save('../../data/Kc167_exp1/local_norm_factor.Kc167_exp1.mat', 'local_norm_factor', 'bad_region_Filter')
