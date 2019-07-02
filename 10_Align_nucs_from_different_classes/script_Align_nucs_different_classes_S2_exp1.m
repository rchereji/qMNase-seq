Figures_folder = 'S2_exp1';
mkdir(Figures_folder)

% Load data
load('../../data/S2_exp1/Fit_results.100-200.NF40.S2_exp1.mat')
[~, max_ids] = arrayfun(@(i) max(NucCount(i,:)), 1:size(NucCount, 1));

% Compute percentage of nucleosomes from each class
edges = 0.5:1:6.5;
Counts = histc(max_ids, edges); Counts = Counts(1:6);
norm_Counts = 100*Counts/sum(Counts);

%%
beforeSite = 1000;
afterSite = 1000;

%%
my_ids = find(max_ids == 1);

myLoc = Loc(my_ids);
myChr = Chr(my_ids);


load('../../data/S2_exp1/local_norm_factor.S2_exp1.mat', 'bad_region_Filter', 'local_norm_factor')
AlignedFilter = align_sites(bad_region_Filter, myChr, myLoc, beforeSite, afterSite);
Filter_score = sum(AlignedFilter, 2);
Filter_score(Filter_score > 0) = NaN;

% Eliminate the genes where the coverage is extremely high or low
badLoci = isnan(Filter_score);

myChr(badLoci) = [];
myLoc(badLoci) = [];


% Align profiles
load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_1min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_01 = nanmean(AlignedProfile);

load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_2min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_02 = nanmean(AlignedProfile);

load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_5min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_05 = nanmean(AlignedProfile);

load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_15min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_15 = nanmean(AlignedProfile);

load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_40min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_40 = nanmean(AlignedProfile);

load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_60min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_60 = nanmean(AlignedProfile);

% Plot average profiles
figure('Position', [50, 50, 450, 250])
hold all
plot([-beforeSite:afterSite]/1000, Avg_counts_01)
plot([-beforeSite:afterSite]/1000, Avg_counts_02)
plot([-beforeSite:afterSite]/1000, Avg_counts_05)
plot([-beforeSite:afterSite]/1000, Avg_counts_15)
plot([-beforeSite:afterSite]/1000, Avg_counts_40)
plot([-beforeSite:afterSite]/1000, Avg_counts_60)
ylim([0, 55])

xlabel('Position relative to peak center (kb)')
ylabel('Nucleosome counts')
title(sprintf('Nucleosomes in class %d: %0.2f%%', 1, norm_Counts(1)))

h = legend({'1 min', '2 min', '5 min', '15 min', '40 min', '60 min'}, 'location', 'EO');
v = get(h,'title');
set(v,'string','Digestion time');
set(gca, 'YGrid', 'on', 'GridLineStyle', '--', 'layer', 'top')
set(gca, 'FontSize', 14)

set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'Nuc_counts_group_1.S2_exp1.eps'));


%%
my_ids = find(max_ids == 2);

myLoc = Loc(my_ids);
myChr = Chr(my_ids);


load('../../data/S2_exp1/local_norm_factor.S2_exp1.mat', 'bad_region_Filter', 'local_norm_factor')
AlignedFilter = align_sites(bad_region_Filter, myChr, myLoc, beforeSite, afterSite);
Filter_score = sum(AlignedFilter, 2);
Filter_score(Filter_score > 0) = NaN;

% Eliminate the genes where the coverage is extremely high or low
badLoci = isnan(Filter_score);

myChr(badLoci) = [];
myLoc(badLoci) = [];


% Align profiles
load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_1min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_01 = nanmean(AlignedProfile);

load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_2min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_02 = nanmean(AlignedProfile);

load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_5min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_05 = nanmean(AlignedProfile);

load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_15min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_15 = nanmean(AlignedProfile);

load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_40min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_40 = nanmean(AlignedProfile);

load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_60min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_60 = nanmean(AlignedProfile);

% Plot average profiles
figure('Position', [50, 50, 450, 250])
hold all
plot([-beforeSite:afterSite]/1000, Avg_counts_01)
plot([-beforeSite:afterSite]/1000, Avg_counts_02)
plot([-beforeSite:afterSite]/1000, Avg_counts_05)
plot([-beforeSite:afterSite]/1000, Avg_counts_15)
plot([-beforeSite:afterSite]/1000, Avg_counts_40)
plot([-beforeSite:afterSite]/1000, Avg_counts_60)
ylim([0, 55])

xlabel('Position relative to peak center (kb)')
ylabel('Nucleosome counts')
title(sprintf('Nucleosomes in class %d: %0.2f%%', 2, norm_Counts(2)))

h = legend({'1 min', '2 min', '5 min', '15 min', '40 min', '60 min'}, 'location', 'EO');
v = get(h,'title');
set(v,'string','Digestion time');
set(gca, 'YGrid', 'on', 'GridLineStyle', '--', 'layer', 'top')
set(gca, 'FontSize', 14)

set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'Nuc_counts_group_2.S2_exp1.eps'));


%%
my_ids = find(max_ids == 3);

myLoc = Loc(my_ids);
myChr = Chr(my_ids);


load('../../data/S2_exp1/local_norm_factor.S2_exp1.mat', 'bad_region_Filter', 'local_norm_factor')
AlignedFilter = align_sites(bad_region_Filter, myChr, myLoc, beforeSite, afterSite);
Filter_score = sum(AlignedFilter, 2);
Filter_score(Filter_score > 0) = NaN;

% Eliminate the genes where the coverage is extremely high or low
badLoci = isnan(Filter_score);

myChr(badLoci) = [];
myLoc(badLoci) = [];


% Align profiles
load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_1min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_01 = nanmean(AlignedProfile);

load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_2min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_02 = nanmean(AlignedProfile);

load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_5min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_05 = nanmean(AlignedProfile);

load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_15min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_15 = nanmean(AlignedProfile);

load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_40min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_40 = nanmean(AlignedProfile);

load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_60min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_60 = nanmean(AlignedProfile);

% Plot average profiles
figure('Position', [50, 50, 450, 250])
hold all
plot([-beforeSite:afterSite]/1000, Avg_counts_01)
plot([-beforeSite:afterSite]/1000, Avg_counts_02)
plot([-beforeSite:afterSite]/1000, Avg_counts_05)
plot([-beforeSite:afterSite]/1000, Avg_counts_15)
plot([-beforeSite:afterSite]/1000, Avg_counts_40)
plot([-beforeSite:afterSite]/1000, Avg_counts_60)
ylim([0, 55])

xlabel('Position relative to peak center (kb)')
ylabel('Nucleosome counts')
title(sprintf('Nucleosomes in class %d: %0.2f%%', 3, norm_Counts(3)))

h = legend({'1 min', '2 min', '5 min', '15 min', '40 min', '60 min'}, 'location', 'EO');
v = get(h,'title');
set(v,'string','Digestion time');
set(gca, 'YGrid', 'on', 'GridLineStyle', '--', 'layer', 'top')
set(gca, 'FontSize', 14)

set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'Nuc_counts_group_3.S2_exp1.eps'));


%%
my_ids = find(max_ids == 4);

myLoc = Loc(my_ids);
myChr = Chr(my_ids);


load('../../data/S2_exp1/local_norm_factor.S2_exp1.mat', 'bad_region_Filter', 'local_norm_factor')
AlignedFilter = align_sites(bad_region_Filter, myChr, myLoc, beforeSite, afterSite);
Filter_score = sum(AlignedFilter, 2);
Filter_score(Filter_score > 0) = NaN;

% Eliminate the genes where the coverage is extremely high or low
badLoci = isnan(Filter_score);

myChr(badLoci) = [];
myLoc(badLoci) = [];


% Align profiles
load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_1min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_01 = nanmean(AlignedProfile);

load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_2min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_02 = nanmean(AlignedProfile);

load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_5min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_05 = nanmean(AlignedProfile);

load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_15min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_15 = nanmean(AlignedProfile);

load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_40min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_40 = nanmean(AlignedProfile);

load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_60min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_60 = nanmean(AlignedProfile);

% Plot average profiles
figure('Position', [50, 50, 450, 250])
hold all
plot([-beforeSite:afterSite]/1000, Avg_counts_01)
plot([-beforeSite:afterSite]/1000, Avg_counts_02)
plot([-beforeSite:afterSite]/1000, Avg_counts_05)
plot([-beforeSite:afterSite]/1000, Avg_counts_15)
plot([-beforeSite:afterSite]/1000, Avg_counts_40)
plot([-beforeSite:afterSite]/1000, Avg_counts_60)
ylim([0, 55])

xlabel('Position relative to peak center (kb)')
ylabel('Nucleosome counts')
title(sprintf('Nucleosomes in class %d: %0.2f%%', 4, norm_Counts(4)))

h = legend({'1 min', '2 min', '5 min', '15 min', '40 min', '60 min'}, 'location', 'EO');
v = get(h,'title');
set(v,'string','Digestion time');
set(gca, 'YGrid', 'on', 'GridLineStyle', '--', 'layer', 'top')
set(gca, 'FontSize', 14)

set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'Nuc_counts_group_4.S2_exp1.eps'));


%%
my_ids = find(max_ids == 5);

myLoc = Loc(my_ids);
myChr = Chr(my_ids);


load('../../data/S2_exp1/local_norm_factor.S2_exp1.mat', 'bad_region_Filter', 'local_norm_factor')
AlignedFilter = align_sites(bad_region_Filter, myChr, myLoc, beforeSite, afterSite);
Filter_score = sum(AlignedFilter, 2);
Filter_score(Filter_score > 0) = NaN;

% Eliminate the genes where the coverage is extremely high or low
badLoci = isnan(Filter_score);

myChr(badLoci) = [];
myLoc(badLoci) = [];


% Align profiles
load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_1min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_01 = nanmean(AlignedProfile);

load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_2min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_02 = nanmean(AlignedProfile);

load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_5min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_05 = nanmean(AlignedProfile);

load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_15min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_15 = nanmean(AlignedProfile);

load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_40min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_40 = nanmean(AlignedProfile);

load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_60min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_60 = nanmean(AlignedProfile);

% Plot average profiles
figure('Position', [50, 50, 450, 250])
hold all
plot([-beforeSite:afterSite]/1000, Avg_counts_01)
plot([-beforeSite:afterSite]/1000, Avg_counts_02)
plot([-beforeSite:afterSite]/1000, Avg_counts_05)
plot([-beforeSite:afterSite]/1000, Avg_counts_15)
plot([-beforeSite:afterSite]/1000, Avg_counts_40)
plot([-beforeSite:afterSite]/1000, Avg_counts_60)
ylim([0, 55])

xlabel('Position relative to peak center (kb)')
ylabel('Nucleosome counts')
title(sprintf('Nucleosomes in class %d: %0.2f%%', 5, norm_Counts(5)))

h = legend({'1 min', '2 min', '5 min', '15 min', '40 min', '60 min'}, 'location', 'EO');
v = get(h,'title');
set(v,'string','Digestion time');
set(gca, 'YGrid', 'on', 'GridLineStyle', '--', 'layer', 'top')
set(gca, 'FontSize', 14)

set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'Nuc_counts_group_5.S2_exp1.eps'));


%%
my_ids = find(max_ids == 6);

myLoc = Loc(my_ids);
myChr = Chr(my_ids);


load('../../data/S2_exp1/local_norm_factor.S2_exp1.mat', 'bad_region_Filter', 'local_norm_factor')
AlignedFilter = align_sites(bad_region_Filter, myChr, myLoc, beforeSite, afterSite);
Filter_score = sum(AlignedFilter, 2);
Filter_score(Filter_score > 0) = NaN;

% Eliminate the genes where the coverage is extremely high or low
badLoci = isnan(Filter_score);

myChr(badLoci) = [];
myLoc(badLoci) = [];


% Align profiles
load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_1min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_01 = nanmean(AlignedProfile);

load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_2min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_02 = nanmean(AlignedProfile);

load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_5min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_05 = nanmean(AlignedProfile);

load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_15min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_15 = nanmean(AlignedProfile);

load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_40min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_40 = nanmean(AlignedProfile);

load('../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_60min.dm6.mat', 'Occ_147');
Norm_occ = cellfun(@(x,y) x./y, Occ_147, local_norm_factor, 'un', 0);
AlignedProfile = align_sites(Norm_occ, myChr, myLoc, beforeSite, afterSite);
Avg_counts_60 = nanmean(AlignedProfile);

% Plot average profiles
figure('Position', [50, 50, 450, 250])
hold all
plot([-beforeSite:afterSite]/1000, Avg_counts_01)
plot([-beforeSite:afterSite]/1000, Avg_counts_02)
plot([-beforeSite:afterSite]/1000, Avg_counts_05)
plot([-beforeSite:afterSite]/1000, Avg_counts_15)
plot([-beforeSite:afterSite]/1000, Avg_counts_40)
plot([-beforeSite:afterSite]/1000, Avg_counts_60)
ylim([0, 55])

xlabel('Position relative to peak center (kb)')
ylabel('Nucleosome counts')
title(sprintf('Nucleosomes in class %d: %0.2f%%', 6, norm_Counts(6)))

h = legend({'1 min', '2 min', '5 min', '15 min', '40 min', '60 min'}, 'location', 'EO');
v = get(h,'title');
set(v,'string','Digestion time');
set(gca, 'YGrid', 'on', 'GridLineStyle', '--', 'layer', 'top')
set(gca, 'FontSize', 14)

set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(Figures_folder, 'Nuc_counts_group_6.S2_exp1.eps'));
