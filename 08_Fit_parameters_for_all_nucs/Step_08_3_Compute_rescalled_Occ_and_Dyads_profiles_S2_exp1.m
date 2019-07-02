%%
load('../../data/S2_exp1/All_Raw_Dyads.100-200.S2_exp1.dm6.mat', 'Dyads')
Occ_147 = extend_dyads_to_fixed_footprint(Dyads, 147);

chrLen = cellfun(@length, Dyads);
% Initialize a filter for the number of nucleosomes (in case of overlapping
% nucs) and a local rescalling factor
filterNoNucs = cell(size(Dyads));
rescallingFactor = cell(size(Dyads));
for c = 1:length(chrLen)
    filterNoNucs{c} = zeros(1, chrLen(c));
    rescallingFactor{c} = nan(1, chrLen(c));
end

load('../../data/S2_exp1/Fit_results.100-200.NF40.S2_exp1.mat', 'Chr', 'Loc', 'O_vector');

noNucs = numel(Chr);
for n = 1 : noNucs   
    filterNoNucs{Chr(n)}(Loc(n)+[-73:73]) = filterNoNucs{Chr(n)}(Loc(n)+[-73:73]) + 1;
    rescallingFactor{Chr(n)}(Loc(n)+[-73:73]) = O_vector(n) / max(Occ_147{Chr(n)}(Loc(n)+[-73:73]));
end

% Set rescallingFactor to NaN where there are overlapping nucleosomes
% (because the correct occupancy at these regions is not clear)
for c = 1:length(chrLen)
    overlappingIdx = filterNoNucs{c} > 1;
    rescallingFactor{c}(overlappingIdx) = nan;
end

% Interpolate the rescallingFactor; create rescallingFactor_v2
filledInd = cell(size(rescallingFactor));
rescallingFactor_v2 = cell(size(rescallingFactor));
for c = 1:length(chrLen)
    [rescallingFactor_v2{c}, filledInd{c}] = fillmissing(rescallingFactor{c}, 'spline','EndValues','nearest');
    disp(c)
end

% Rescale the Dyads and Occ profiles according to the local rescalling
% factor
rescalledDyads = cell(size(Dyads));
rescalledOcc = cell(size(Dyads));
for c = 1:length(chrLen)
    rescalledDyads{c} = rescallingFactor_v2{c} .* Dyads{c};
    rescalledOcc{c} = rescallingFactor_v2{c} .* Occ_147{c};
end

% Fix a problem due to spline interpolation: some regions got a high
% rescalling factor and consequently an Occ > 1
% Mark these bad regions
badRegion = cell(size(rescalledOcc));
for c = 1:length(chrLen)
    badRegion{c} = (rescalledOcc{c} > 1);
    rescalledOcc{c}(badRegion{c}) = 1;
end

% Save the rescalled profiles
save('../../data/S2_exp1/Rescalled_Dyads_and_Occ.100-200.S2_exp1.mat', ...
    'rescalledDyads', 'rescalledOcc', 'badRegion');
