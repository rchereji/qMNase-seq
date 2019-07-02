sampleLabel = {'1min', '2min', '5min', '15min', '40min', '60min'};
noSamples = numel(sampleLabel);

%% S2_exp1 reads
pathName = '../../data/S2_exp1';

s = 1;
load(fullfile(pathName, ['Raw_Dyads.100-200.S2_exp1_', sampleLabel{s}, '.dm6.mat']), 'Dyads', 'noReadsPerChr', 'totalNoReads', 'chrName', 'chrLen')
Dyads_all_datasets = Dyads;
noReadsPerChr_all_datasets = noReadsPerChr;
totalNoReads_all_datasets = totalNoReads;

for s = 2:noSamples
    load(fullfile(pathName, ['Raw_Dyads.100-200.S2_exp1_', sampleLabel{s}, '.dm6.mat']), 'Dyads', 'noReadsPerChr', 'totalNoReads')
    Dyads_all_datasets = cellfun(@(x, y) x+y, Dyads, Dyads_all_datasets, 'un', 0);
    noReadsPerChr_all_datasets = noReadsPerChr_all_datasets + noReadsPerChr;
    totalNoReads_all_datasets = totalNoReads_all_datasets + totalNoReads;
end

Dyads = Dyads_all_datasets;
totalNoReads = totalNoReads_all_datasets;
noReadsPerChr = noReadsPerChr_all_datasets;
save(fullfile(pathName, 'All_Raw_Dyads.100-200.S2_exp1.dm6.mat'), 'Dyads', 'noReadsPerChr', 'totalNoReads', 'chrName', 'chrLen')

%% S2_exp2 reads
pathName = '../../data/S2_exp2';

s = 1;
load(fullfile(pathName, ['Raw_Dyads.100-200.S2_exp2_seq1_', sampleLabel{s}, '.dm6.mat']), 'Dyads', 'noReadsPerChr', 'totalNoReads', 'chrName', 'chrLen')
Dyads_all_datasets = Dyads;
noReadsPerChr_all_datasets = noReadsPerChr;
totalNoReads_all_datasets = totalNoReads;

for s = 2:noSamples
    load(fullfile(pathName, ['Raw_Dyads.100-200.S2_exp2_seq1_', sampleLabel{s}, '.dm6.mat']), 'Dyads', 'noReadsPerChr', 'totalNoReads')
    Dyads_all_datasets = cellfun(@(x, y) x+y, Dyads, Dyads_all_datasets, 'un', 0);
    noReadsPerChr_all_datasets = noReadsPerChr_all_datasets + noReadsPerChr;
    totalNoReads_all_datasets = totalNoReads_all_datasets + totalNoReads;
end

Dyads = Dyads_all_datasets;
totalNoReads = totalNoReads_all_datasets;
noReadsPerChr = noReadsPerChr_all_datasets;
save(fullfile(pathName, 'All_Raw_Dyads.100-200.S2_exp2.dm6.mat'), 'Dyads', 'noReadsPerChr', 'totalNoReads', 'chrName', 'chrLen')

%% Kc167_exp1 reads
pathName = '../../data/Kc167_exp1';

s = 1;
load(fullfile(pathName, ['Raw_Dyads.100-200.Kc167_exp1_seq1_', sampleLabel{s}, '.dm6.mat']), 'Dyads', 'noReadsPerChr', 'totalNoReads', 'chrName', 'chrLen')
Dyads_all_datasets = Dyads;
noReadsPerChr_all_datasets = noReadsPerChr;
totalNoReads_all_datasets = totalNoReads;

for s = 2:noSamples
    load(fullfile(pathName, ['Raw_Dyads.100-200.Kc167_exp1_seq1_', sampleLabel{s}, '.dm6.mat']), 'Dyads', 'noReadsPerChr', 'totalNoReads')
    Dyads_all_datasets = cellfun(@(x, y) x+y, Dyads, Dyads_all_datasets, 'un', 0);
    noReadsPerChr_all_datasets = noReadsPerChr_all_datasets + noReadsPerChr;
    totalNoReads_all_datasets = totalNoReads_all_datasets + totalNoReads;
end

Dyads = Dyads_all_datasets;
totalNoReads = totalNoReads_all_datasets;
noReadsPerChr = noReadsPerChr_all_datasets;
save(fullfile(pathName, 'All_Raw_Dyads.100-200.Kc167_exp1.dm6.mat'), 'Dyads', 'noReadsPerChr', 'totalNoReads', 'chrName', 'chrLen')
