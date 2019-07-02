function count_Scer_spikeins(bamFilename, Lmin, Lmax)
% This function loads the alignments from a BAM file containing the spikeins, 
% and computes the number of spikeins with the lengths between Lmin and Lmax 
%
% Inputs:
% bamFilename  - BAM file containing the genomic alignments of the
%                sequenced DNA fragments
% Lmin / Lmax  - DNA sizes to be considered (e.g. 120-160 bp)
% 
% Output: 
% A MATLAB file with the no. of reads detected on each chromosome, with the
% length between Lmin and Lmax
%
% Example:
% count_Scer_spikeins('myData.Scer_spikeins.bam', 100, 200);

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

% Get chrom names and lengths
InfoStruct = baminfo(bamFilename);
chrName = {InfoStruct.SequenceDictionary.SequenceName};
chrLen = [InfoStruct.SequenceDictionary.SequenceLength];
noChr = numel(chrName);

%% Process all chromosomes
noSpikeinsPerChr = zeros(1, noChr);
totalNoSpikeins = 0;

for chr = 1 : noChr
    % Create BioMap object
    bm = BioMap(bamFilename, 'SelectReference', chrName{chr});
    
    % Eliminate the reads with low MAPping Quality (MAPQ). 
    % MAPQ = ?10 log10 Prob(mapping position is wrong)
    bm_filtered = getSubset(bm, getMappingQuality(bm) >= 30);   % Use only the reads with QMAP >= 30
    
    % Filter out the reads which fall outside the reference, if any alignment artifacts are detected
    bm_filtered = getSubset(bm_filtered, getStop(bm_filtered) <= chrLen(chr));
    
    Indices = filterByFlag(bm_filtered, 'pairedInMap', true, 'strandQuery', 0);
    bm_filtered_Watson = getSubset(bm_filtered, Indices);
    
    Indices = filterByFlag(bm_filtered, 'pairedInMap', true, 'strandQuery', 1);
    bm_filtered_Crick = getSubset(bm_filtered, Indices);
    
	[~, Watson_Idx, Crick_Idx] = intersect(bm_filtered_Watson.Header, bm_filtered_Crick.Header);
    
    % Compute fragment length histogram
    leftBP = getStart(bm_filtered_Watson, Watson_Idx);
    rightBP = getStop(bm_filtered_Crick, Crick_Idx);
    
    fragmentLengths = rightBP - leftBP + 1;
    
    goodLengthFilter = ((fragmentLengths >= Lmin) & (fragmentLengths <= Lmax));
    noSpikeinsPerChr(chr) = sum(goodLengthFilter);
    totalNoSpikeins = totalNoSpikeins + noSpikeinsPerChr(chr);
end

% Save the numbers of detected spike-ins
save(fullfile(pathName, sprintf('Scer_spikeins.%d-%d.%s.mat', Lmin, Lmax, erase(setName, '.Scer_spikeins'))), ...
    'totalNoSpikeins', 'noSpikeinsPerChr', 'chrName', 'chrLen');
