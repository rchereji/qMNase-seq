% S2_exp1
pathName = '../../data/S2_exp1';
allFiles = dir([pathName, '/*.bam']);

for f = 1:length(allFiles)
    bamFilename = fullfile(pathName, allFiles(f).name);
    get_DNA_fragment_lengths(bamFilename);
end

% S2_exp2
pathName = '../../data/S2_exp2';
allFiles = dir([pathName, '/*.bam']);

for f = 1:length(allFiles)
    bamFilename = fullfile(pathName, allFiles(f).name);
    get_DNA_fragment_lengths(bamFilename);
end

% Kc167_exp1
pathName = '../../data/Kc167_exp1';
allFiles = dir([pathName, '/*.bam']);

for f = 1:length(allFiles)
    bamFilename = fullfile(pathName, allFiles(f).name);
    get_DNA_fragment_lengths(bamFilename);
end
