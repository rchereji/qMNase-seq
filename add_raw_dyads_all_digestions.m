sampleLabel = {'1min', '2min', '5min', '15min', '40min', '60min'};
noSamples = numel(sampleLabel);

s = 1;
load(['data/Raw_Dyads.100-200.S2_', sampleLabel{s}, '.dm6.mat'], 'Dyads', 'noReadsPerChr', 'totalNoReads', 'chrName', 'chrLen')
Dyads_all_datasets = Dyads;
noReadsPerChr_all_datasets = noReadsPerChr;
totalNoReads_all_datasets = totalNoReads;

for s = 2:noSamples
    load(['data/Raw_Dyads.100-200.S2_', sampleLabel{s}, '.dm6.mat'], 'Dyads', 'noReadsPerChr', 'totalNoReads')
    Dyads_all_datasets = cellfun(@(x, y) x+y, Dyads, Dyads_all_datasets, 'un', 0);
    noReadsPerChr_all_datasets = noReadsPerChr_all_datasets + noReadsPerChr;
    totalNoReads_all_datasets = totalNoReads_all_datasets + totalNoReads;
end

%%
Dyads = Dyads_all_datasets;
totalNoReads = totalNoReads_all_datasets;
noReadsPerChr = noReadsPerChr_all_datasets;
save('data/All_Dyads.100-200.S2_cells.dm6.mat', 'Dyads', 'noReadsPerChr', 'totalNoReads', 'chrName', 'chrLen')
