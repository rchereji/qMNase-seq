% Size selection: use reads of 100-200 bp to count the nucleosomal fragments
Lmin = 100;
Lmax = 200;

% Read the available data files (BAM format)
% S2_exp1
pathName = '../../data/S2_exp1';
allFiles = dir([pathName, '/*.bam']);

for f = 1:length(allFiles)
    bamFilename = fullfile(pathName, allFiles(f).name);
    get_raw_dyads_and_occupancy(bamFilename, Lmin, Lmax);
end


% S2_exp2
pathName = '../../data/S2_exp2';
allFiles = dir([pathName, '/*.bam']);

for f = 1:length(allFiles)
    bamFilename = fullfile(pathName, allFiles(f).name);
    get_raw_dyads_and_occupancy(bamFilename, Lmin, Lmax);
end


% Kc167 reads
pathName = '../../data/Kc167_exp1';
allFiles = dir([pathName, '/*.bam']);

for f = 1:length(allFiles)
    bamFilename = fullfile(pathName, allFiles(f).name);
    get_raw_dyads_and_occupancy(bamFilename, Lmin, Lmax);
end
