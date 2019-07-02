function alignedProfiles = align_sites(Profile, Chr, Loc, beforeSite, afterSite)

% Compute the chromosome sizes
chrLen = cellfun(@length, Profile);

noSites = length(Chr);
alignedProfiles = nan(noSites, 1 + beforeSite + afterSite);
for s = 1:noSites
    if ~isnan(Loc(s))
        leftEdge = max([Loc(s) - beforeSite, 1]);
        rightEdge = min([Loc(s) + afterSite, chrLen(Chr(s))]);
        alignedProfiles(s, beforeSite + 1 - (Loc(s) - leftEdge)...
            : beforeSite + 1 + (rightEdge - Loc(s))) = ...
            Profile{Chr(s)}(leftEdge : rightEdge);
    end
end
