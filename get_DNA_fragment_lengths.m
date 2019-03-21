function [frag_length, frag_count] = get_DNA_fragment_lengths(bamFilename, noReadsToStop)
% This function computes the read length distribution from the BAM file
%
% Inputs:
% bamFilename   - BAM file containing the genomic alignments of the
%                 sequenced DNA fragments
% noReadsToStop - the minimum number of reads to be used to generate the 
%                 histogram (default value: 1,000,000)
%
% Outputs:
% frag_length - length of a DNA fragment, between 1 and 1000 bp
% frag_count  - the number of reads corresponding to each DNA size
%
% Example:
% [frag_length, frag_count] = get_DNA_fragment_lengths('myData.bam');

warning('off', 'bioinfo:saminfo:InvalidTagField')

% Check the input, get setName
if nargin == 0
    error('You didn''t provide the data file!')
else
    if (exist(bamFilename, 'file') ~= 2)
        error('File "%s" does not exist in the current folder!', bamFilename)
    end
    
    [~,setName,ext] = fileparts(bamFilename);
    if ~strcmpi(ext, '.bam')
        error('File "%s" is not a BAM file!', bamFilename)
    end
end

if exist('noReadsToStop', 'var')
    % make sure the parameter is numeric, not a string
    if ischar(noReadsToStop)
        noReadsToStop = str2double(noReadsToStop);
    end
else
    % if noReadsToStop does not exist, set default value
    noReadsToStop = 1000000;
end

% Get chrom names and lengths
InfoStruct = baminfo(bamFilename);
chrName = {InfoStruct.SequenceDictionary.SequenceName};
chrLen = [InfoStruct.SequenceDictionary.SequenceLength];

% Process all chromosomes
noChr = numel(chrName);
all_fragment_lengths = [];

for chr = 1 : noChr
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
    
    all_fragment_lengths = [all_fragment_lengths; fragmentLengths];
    
    % Stop when at least 100,000 proper paired-end reads have been analyzed
    % If all the fragments need to be considered, just comment out the following IF statement
    if sum(all_fragment_lengths) > noReadsToStop
        break
    end
end

frag_length = [1:1000]';
frag_count = arrayfun(@(z) sum(ismember(all_fragment_lengths, z)), frag_length);

save(sprintf('Length_histogram.%s.mat', setName), 'frag_length', 'frag_count');

% Plot the figure
figure('Position', [50, 50, 400, 300]);
plot(frag_length, 100*frag_count/sum(frag_count), 'LineWidth', 1);
xlim([0 500])
set(gca, 'XTick', 0:100:500, 'XGrid', 'on', 'GridLineStyle', '--')
set(gca, 'XMinorTick', 'on', 'TickLength', [0.03 0.025])
set(gca, 'FontSize', 12)
xlabel('Fragment length (bp)', 'FontSize', 14)
ylabel('Percentage (%)', 'FontSize', 14)
title({'Length histogram'; sprintf('Sample: %s', setName)}, 'interpreter', 'none', 'FontSize', 14)
box off

% Create inset with mono-nuc data
axes('Position',[.49 .42 .38 .4])
plot(frag_length, 100*frag_count/sum(frag_count), 'LineWidth', 2);
xlim([100 200])
set(gca, 'XTick', 100:20:200, 'XGrid', 'on', 'GridLineStyle', '--')
set(gca, 'FontSize', 10)
% title({'Zoom-in view:'; '100 bp - 200 bp'}, 'FontSize', 12)

% Save the figure as an EPS file
set(gcf, 'PaperPositionMode', 'auto')
print(gcf, sprintf('Length_histogram.%s.eps', setName), '-depsc', '-painters');
