function [Dyads, Occ] = get_raw_dyads_and_occupancy(bamFilename, Lmin, Lmax)
% This function loads all the alignments from a BAM file, and computes the
% coverage profile (Occ) and the distribution of the fragment centers (Dyads)
%
% Inputs:
% bamFilename - BAM file containing the genomic alignments of the
% sequenced DNA fragments
% Lmin / Lmax  - DNA sizes to be considered (e.g. 120-160 bp)
% 
% Output: 
% Dyads - cell array (one element for each chromosome) containing the
%         number of nucleosome centers (dyads) that were detected at each
%         genomic position 
% Occ   - cell array containing the number of DNA fragments that cover
%         each genomic position (i.e. occupancy / coverage profile)
%
% Example:
% [Dyads, Occ] = get_dyads_and_occupancy('myData.bam', 120, 160);

warning('off', 'bioinfo:saminfo:InvalidTagField')

% Check the input, get setName
if nargin == 0
    error('You didn''t provide the data file!')
else
    if (exist(bamFilename, 'file') ~= 2)
        error('File "%s" does not exist in the current folder!', bamFilename)
    end
    
    [pathName, setName, ext] = fileparts(bamFilename);
    if ~strcmpi(ext, '.bam')
        error('File "%s" is not a BAM file!', bamFilename)
    end
    
    % Check Lmin argument
    if exist('Lmin', 'var')
        % make sure the parameter is numeric, not a string
        if ischar(Lmin)
            Lmin = str2double(Lmin);
        end
    else
        % if Lmin does not exist, set default value
        Lmin = 100;
    end

    % Check Lmax argument
    if exist('Lmax', 'var')
        % make sure the parameter is numeric, not a string
        if ischar(Lmax)
            Lmax = str2double(Lmax);
        end
    else
        % if Lmax does not exist, set default value
        Lmax = 200;
    end
end

tic

% Get chrom names and lengths
InfoStruct = baminfo(bamFilename);
chrName = {InfoStruct.SequenceDictionary.SequenceName};
chrLen = [InfoStruct.SequenceDictionary.SequenceLength];
noChr = numel(chrName);

% Initialize the Occ & Dyads cell arrays
Occ = cell(1, noChr);
Dyads = cell(1, noChr);

fprintf('Starting to import data from "%s".\n', bamFilename)
fprintf('Size selection: %d <= L <= %d\n', Lmin, Lmax)

% Process all chromosomes; parallel for loop, using the maximum number of workers/threads
noReadsPerChr = zeros(1, noChr);

for chr = 1 : noChr
% parfor chr = 1 : noChr   % Multi-threaded version for speed improvement

    % Create BioMap object
    bm = BioMap(bamFilename, 'SelectReference', chrName{chr});
    
    % Eliminate the reads with very low MAPping Quality (MAPQ = 0). 
    % MAPQ = ?10 log10 Prob(mapping position is wrong)
    bm = getSubset(bm, getMappingQuality(bm) > 0);
    
    % Filter out the reads which fall outside the reference, if any alignment artifacts are detected
    bm = getSubset(bm, getStop(bm) <= chrLen(chr));
    
    % Get the reads that map to the Watson strand
    Indices = filterByFlag(bm, 'pairedInMap', true, 'strandQuery', 0);
    bm_filtered_Watson = getSubset(bm, Indices);
    
    % Get the reads that map to the Crick strand
    Indices = filterByFlag(bm, 'pairedInMap', true, 'strandQuery', 1);
    bm_filtered_Crick = getSubset(bm, Indices);
    
    % Match the paired-end reads (1 read mapped to Watson strand and 1 read mapped to Crick strand)
	[~, Watson_Idx, Crick_Idx] = intersect(bm_filtered_Watson.Header, bm_filtered_Crick.Header);
    
    % Compute the fragment lengths
    leftBP = getStart(bm_filtered_Watson, Watson_Idx);
    rightBP = getStop(bm_filtered_Crick, Crick_Idx);
    fragmentLengths = rightBP - leftBP + 1;
    
    % Size selection, acording to the specified limits (Lmin, Lmax)
    goodInd = ((fragmentLengths >= Lmin) & (fragmentLengths <= Lmax));
    leftBP = leftBP(goodInd);
    rightBP = rightBP(goodInd);
    
    noReadsPerChr(chr) = numel(leftBP);
    
    % Construct the distribution of fragment centers (Dyads)
    Centers = round((rightBP + leftBP)/2);
    uniqueCenters = unique(Centers);
    [~, Index] = ismember(Centers, uniqueCenters);
    NumberUniqueCenter = histc(Index, 1:numel(uniqueCenters));
    
    Dyads{chr} = zeros(1, chrLen(chr));
    Dyads{chr}(uniqueCenters) = NumberUniqueCenter;
    
    % Compute the coverage/occupancy (Occ distribution)
    uniqueLeftBP = unique(leftBP);
    [~, Index] = ismember(leftBP, uniqueLeftBP);
    NumberUniqueleftBP = histc(Index, 1:numel(uniqueLeftBP));
    
    OccDerivative = zeros(1, chrLen(chr) + 1); % add 1 position for the case when some reads have the right end exactly at the end of chr
    OccDerivative(uniqueLeftBP) = NumberUniqueleftBP;
    
    uniqueRightBP = unique(rightBP);
    [~, Index] = ismember(rightBP, uniqueRightBP);
    NumberUniqueRightBP = histc(Index, 1:numel(uniqueRightBP));
    
    OccDerivative(uniqueRightBP + 1) = OccDerivative(uniqueRightBP + 1) - NumberUniqueRightBP';
    tmp = cumsum(OccDerivative);
    
    Occ{chr} = tmp(1 : end-1);
    
    fprintf('Chromosome %s done.\n', chrName{chr});
end

% Create a log string
infoStr = sprintf('Sequencing depth report\nDataset: %s\nSize selection: %d-%d bp\n\nReference    Number of Pairs    Density (pairs/bp)\n', bamFilename, Lmin, Lmax);

for chr = 1 : noChr
    % Fill in the log string
    infoStr = sprintf('%s%8s\t\t%d\t\t\t\t%0.6f \n', infoStr,...
        chrName{chr}, noReadsPerChr(chr), double(noReadsPerChr(chr) / chrLen(chr)));
end 
totalNoReads = sum(noReadsPerChr);
infoStr = sprintf('%s\n%8s\t\t%d\t\t\t\t%0.6f \n', infoStr,...
    'Total', totalNoReads, double(totalNoReads / sum(chrLen)));

% Save the Occupancy and Dyads cell arrays
save(fullfile(pathName, sprintf('Raw_Occ.%d-%d.%s.mat', Lmin, Lmax, setName)), 'Occ', 'totalNoReads', 'noReadsPerChr', 'chrName', 'chrLen');
save(fullfile(pathName, sprintf('Raw_Dyads.%d-%d.%s.mat', Lmin, Lmax, setName)), 'Dyads', 'totalNoReads', 'noReadsPerChr', 'chrName', 'chrLen');

% Save the log file with the sequencing depth
fileID = fopen(fullfile(pathName, sprintf('Sequencing_depth_log.%d-%d.%s.txt', Lmin, Lmax, setName)), 'w');
fprintf(fileID, '%s', infoStr);
fclose(fileID);

tStop = toc;
if tStop > 3600
    fprintf('File "%s.bam" was successfully processed in %d h, %d min, and %0.0f s.\n\n', setName, floor(tStop/3600), floor(mod(tStop, 3600)/60), mod(tStop, 60))
else
    fprintf('File "%s.bam" was successfully processed in %d min and %0.0f s.\n\n', setName, floor(tStop/60), mod(tStop, 60))
end
