function Occ = extend_dyads_to_fixed_footprint(Dyads, L)
% Converts a cell array containing the nucleosome dyad distribution into a
% cell array containing the coverage profile obtained by stacking extended
% footprints, symmetricaly extended from the dyad position (nucleosome
% occupancy profiles, if L = 147 bp).
%
% Inputs:
%   Dyads - cell array, each cell containing the nucleosome dyad counts
%           for an individual chromosome
%   L     - footprint of the particle (typical size for nucleosomes = 147)
%
% Output:
%   Occ   - cell array, each cell containing the nucleosome occupancy
%           for an individual chromosome
%
% Examples:
% Occ = extend_dyads_to_fixed_footprint(Dyads, 147) % compute nucleosome occupancy
% Occ = extend_dyads_to_fixed_footprint(Dyads, 101) % extend dyads to 101 bp 
% (to emphasize the linkers and better distinguish neighboring nucleosomes 
% - useful for visualization purposes, and for detection of the typical 
% positions of nucleosomes)

halfFootprintSize = floor(L/2);
noChr = numel(Dyads);
Occ = cell(size(Dyads));
for c = 1:noChr
    Occ{c} = filter(ones(1,L), 1, [Dyads{c}(halfFootprintSize + 1 : end), zeros(1, halfFootprintSize)]);
end
