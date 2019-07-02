clear
load('AT_GC_content.mat', 'GC_content')
windowWidth = 201;
AvgGCcontent = cellfun(@(x) movmean(x, windowWidth, 'omitnan'), ...
    GC_content, 'un', 0); 

% Create bins og GC-content
GCbinCenters = 0:0.02:1;
GCbinEdges = -0.01:0.02:1.01;

% The AT-content corresponding to the same bins
ATbinCenters = 1 - GCbinCenters;

% Length of the chromosomes
chrLen = cellfun(@length, GC_content);
noChr = numel(chrLen);

%% Compute GC-content, AT-content
%% S2_exp1
% Identify the Dyads files and compute the GC-content
pathName = '../../data/S2_exp1';
dyadDir = dir([pathName, '/Raw_Dyads.*.mat']);

noFiles = numel(dyadDir);
for f = 1:noFiles
    dyadsFilename = dyadDir(f,1).name;
    Dataset = dyadsFilename(11:end-4);
    
    load(fullfile(pathName, dyadsFilename), 'Dyads')
    binCounts = zeros(size(GCbinCenters));
    
    for chr = 1:noChr
        locPerChr = find(Dyads{chr} ~= 0);
        no_dyads = Dyads{chr}(locPerChr);
        badRegion = (no_dyads > quantile(no_dyads, 0.99));
        locPerChr(badRegion) = [];
        no_dyads(badRegion) = [];
        
        noPeaks = numel(locPerChr);
        
        for k = 1:noPeaks
            gc = AvgGCcontent{chr}(locPerChr(k));
            bin = find(GCbinEdges > gc, 1, 'first') - 1;
            binCounts(bin) = binCounts(bin) + no_dyads(k);
        end
        disp(chr)
    end
    
    save(fullfile(pathName, ['AT_GC_distribution.', Dataset,'.mat']), ...
        'ATbinCenters', 'GCbinCenters', 'binCounts', 'windowWidth');
    fprintf('File %s done.\n', dyadsFilename)
end


%% S2_exp2
% Identify the Dyads files and compute the GC-content
pathName = '../../data/S2_exp2';
dyadDir = dir([pathName, '/Raw_Dyads.*.mat']);

noFiles = numel(dyadDir);
for f = 1:noFiles
    dyadsFilename = dyadDir(f,1).name;
    Dataset = dyadsFilename(11:end-4);
    
    load(fullfile(pathName, dyadsFilename), 'Dyads')
    binCounts = zeros(size(GCbinCenters));
    
    for chr = 1:noChr
        locPerChr = find(Dyads{chr} ~= 0);
        no_dyads = Dyads{chr}(locPerChr);
        badRegion = (no_dyads > quantile(no_dyads, 0.99));
        locPerChr(badRegion) = [];
        no_dyads(badRegion) = [];
        
        noPeaks = numel(locPerChr);
        
        for k = 1:noPeaks
            gc = AvgGCcontent{chr}(locPerChr(k));
            bin = find(GCbinEdges > gc, 1, 'first') - 1;
            binCounts(bin) = binCounts(bin) + no_dyads(k);
        end
        disp(chr)
    end
    
    save(fullfile(pathName, ['AT_GC_distribution.', Dataset,'.mat']), ...
        'ATbinCenters', 'GCbinCenters', 'binCounts', 'windowWidth');
    fprintf('File %s done.\n', dyadsFilename)
end


%% Kc167_exp1
% Identify the Dyads files and compute the GC-content
pathName = '../../data/Kc167_exp1';
dyadDir = dir([pathName, '/Raw_Dyads.*.mat']);

noFiles = numel(dyadDir);
for f = 1:noFiles
    dyadsFilename = dyadDir(f,1).name;
    Dataset = dyadsFilename(11:end-4);
    
    load(fullfile(pathName, dyadsFilename), 'Dyads')
    binCounts = zeros(size(GCbinCenters));
    
    for chr = 1:noChr
        locPerChr = find(Dyads{chr} ~= 0);
        no_dyads = Dyads{chr}(locPerChr);
        badRegion = (no_dyads > quantile(no_dyads, 0.99));
        locPerChr(badRegion) = [];
        no_dyads(badRegion) = [];
        
        noPeaks = numel(locPerChr);
        
        for k = 1:noPeaks
            gc = AvgGCcontent{chr}(locPerChr(k));
            bin = find(GCbinEdges > gc, 1, 'first') - 1;
            binCounts(bin) = binCounts(bin) + no_dyads(k);
        end
        disp(chr)
    end
    
    save(fullfile(pathName, ['AT_GC_distribution.', Dataset,'.mat']), ...
        'ATbinCenters', 'GCbinCenters', 'binCounts', 'windowWidth');
    fprintf('File %s done.\n', dyadsFilename)
end

