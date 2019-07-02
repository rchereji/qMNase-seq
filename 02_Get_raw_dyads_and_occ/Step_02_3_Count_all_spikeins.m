% Size selection: use reads of 100-200 bp to count the nucleosomal fragments
Lmin = 100;
Lmax = 200;

% Count spike-ins
% S2_exp1
pathName = '../../data/S2_exp1/spike-ins';
allFiles = dir([pathName, '/*.bam']);

for f = 1:length(allFiles)
    bamFilename = fullfile(pathName, allFiles(f).name);
    count_Scer_spikeins(bamFilename, Lmin, Lmax);
end

% S2_exp2
pathName = '../../data/S2_exp2/spike-ins';
allFiles = dir([pathName, '/*.bam']);

for f = 1:length(allFiles)
    bamFilename = fullfile(pathName, allFiles(f).name);
    count_Scer_spikeins(bamFilename, Lmin, Lmax);
end

% Kc167_exp1
pathName = '../../data/Kc167_exp1/spike-ins';
allFiles = dir([pathName, '/*.bam']);

for f = 1:length(allFiles)
    bamFilename = fullfile(pathName, allFiles(f).name);
    count_Scer_spikeins(bamFilename, Lmin, Lmax);
end
