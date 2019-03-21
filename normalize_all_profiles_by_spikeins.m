%% Load Spike-ins
sampleLabel = {'1min', '2min', '5min', '15min', '40min', '60min'};
noSamples = numel(sampleLabel);

normFactor = nan(1, noSamples);
for s = 1:noSamples
    load(['data/Scer_spikeins.100-200.S2_', sampleLabel{s}, '.mat'], 'totalNoSpikeins')
	normFactor(s) = 10000/totalNoSpikeins;
end

%% Normalize profiles by spike-in counts
for s = 1:noSamples
    load(['data/Raw_Dyads.100-200.S2_', sampleLabel{s}, '.dm6.mat'], 'Dyads', 'totalNoReads', 'noReadsPerChr', 'chrName', 'chrLen')
    totalNoReads = totalNoReads * normFactor(s);
    noReadsPerChr = noReadsPerChr * normFactor(s);
    Dyads = cellfun(@(x) x*normFactor(s), Dyads, 'un', 0);
    Occ_101 = extend_dyads_to_fixed_footprint(Dyads, 101);
    Occ_147 = extend_dyads_to_fixed_footprint(Dyads, 147);
    
    save(['data/Spikein_normalized_profiles.100-200.S2_', sampleLabel{s}, '.dm6.mat'], 'Dyads', 'Occ_101', 'Occ_147', 'totalNoReads', 'noReadsPerChr', 'chrName', 'chrLen')
    disp(s)
end
