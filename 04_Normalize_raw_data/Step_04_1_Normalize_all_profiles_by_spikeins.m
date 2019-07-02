%% Load Spike-ins
sampleLabel = {'1min', '2min', '5min', '15min', '40min', '60min'};
noSamples = numel(sampleLabel);

%% Normalize S2 profiles
% S2_exp1
normFactor = nan(1, noSamples);
for s = 1:noSamples
    load(['../../data/S2_exp1/spike-ins/Scer_spikeins.100-200.S2_exp1_', sampleLabel{s}, '.mat'], 'totalNoSpikeins')
	normFactor(s) = 10000/totalNoSpikeins;
end

for s = 1:noSamples
    load(['../../data/S2_exp1/Raw_Dyads.100-200.S2_exp1_', sampleLabel{s}, '.dm6.mat'], 'Dyads', 'totalNoReads', 'noReadsPerChr', 'chrName', 'chrLen')
    totalNoReads = totalNoReads * normFactor(s);
    noReadsPerChr = noReadsPerChr * normFactor(s);
    Dyads = cellfun(@(x) x*normFactor(s), Dyads, 'un', 0);
    Occ_101 = extend_dyads_to_fixed_footprint(Dyads, 101);
    Occ_147 = extend_dyads_to_fixed_footprint(Dyads, 147);
    
    save(['../../data/S2_exp1/Spikein_normalized_profiles.100-200.S2_exp1_', sampleLabel{s}, '.dm6.mat'], 'Dyads', 'Occ_101', 'Occ_147', 'totalNoReads', 'noReadsPerChr', 'chrName', 'chrLen')
    disp(s)
end


% S2_exp2
normFactor = nan(1, noSamples);
for s = 1:noSamples
    load(['../../data/S2_exp2/spike-ins/Scer_spikeins.100-200.S2_exp2_seq1_', sampleLabel{s}, '.mat'], 'totalNoSpikeins')
	normFactor(s) = 10000/totalNoSpikeins;
end

for s = 1:noSamples
    load(['../../data/S2_exp2/Raw_Dyads.100-200.S2_exp2_seq1_', sampleLabel{s}, '.dm6.mat'], 'Dyads', 'totalNoReads', 'noReadsPerChr', 'chrName', 'chrLen')
    totalNoReads = totalNoReads * normFactor(s);
    noReadsPerChr = noReadsPerChr * normFactor(s);
    Dyads = cellfun(@(x) x*normFactor(s), Dyads, 'un', 0);
    Occ_101 = extend_dyads_to_fixed_footprint(Dyads, 101);
    Occ_147 = extend_dyads_to_fixed_footprint(Dyads, 147);
    
    save(['../../data/S2_exp2/Spikein_normalized_profiles.100-200.S2_exp2_seq1_', sampleLabel{s}, '.dm6.mat'], 'Dyads', 'Occ_101', 'Occ_147', 'totalNoReads', 'noReadsPerChr', 'chrName', 'chrLen')
    disp(s)
end


%% Normalize Kc167 profiles
% Kc167_exp1
normFactor = nan(1, noSamples);
for s = 1:noSamples
    load(['../../data/Kc167_exp1/spike-ins/Scer_spikeins.100-200.Kc167_exp1_seq1_', sampleLabel{s}, '.mat'], 'totalNoSpikeins')
	normFactor(s) = 100000/totalNoSpikeins;
end

for s = 1:noSamples
    load(['../../data/Kc167_exp1/Raw_Dyads.100-200.Kc167_exp1_seq1_', sampleLabel{s}, '.dm6.mat'], 'Dyads', 'totalNoReads', 'noReadsPerChr', 'chrName', 'chrLen')
    totalNoReads = totalNoReads * normFactor(s);
    noReadsPerChr = noReadsPerChr * normFactor(s);
    Dyads = cellfun(@(x) x*normFactor(s), Dyads, 'un', 0);
    Occ_101 = extend_dyads_to_fixed_footprint(Dyads, 101);
    Occ_147 = extend_dyads_to_fixed_footprint(Dyads, 147);
    
    save(['../../data/Kc167_exp1/Spikein_normalized_profiles.100-200.Kc167_exp1_seq1_', sampleLabel{s}, '.dm6.mat'], 'Dyads', 'Occ_101', 'Occ_147', 'totalNoReads', 'noReadsPerChr', 'chrName', 'chrLen')
    disp(s)
end