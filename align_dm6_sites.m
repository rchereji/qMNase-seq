function AlignedProfiles = align_dm6_sites(Profile, Chr, Loc, beforeSite, afterSite)

% Load the chromosome sizes
load('data/chrLen_dm6.mat', 'chrLen')

noSites = length(Chr);
AlignedProfiles = nan(noSites, 1 + beforeSite + afterSite);
for s = 1:noSites
    if ~isnan(Loc(s))
        leftEdge = max([Loc(s) - beforeSite, 1]);
        rightEdge = min([Loc(s) + afterSite, chrLen(Chr(s))]);
        AlignedProfiles(s, beforeSite + 1 - (Loc(s) - leftEdge)...
            : beforeSite + 1 + (rightEdge - Loc(s))) = ...
            Profile{Chr(s)}(leftEdge : rightEdge);
    end
end
