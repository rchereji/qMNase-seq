function write_profile_to_WIG_file(cellArray, chrName, WIG_filename)
% Write a cell array into a WIG file
%
% Inputs:
%   cellArray    - cell array, each cell containing the a profile
%                  corresponding to a chromosome
%   chrName      - chromosome name
%   WIG_filename - name of the WIG file where to save the cell array
%
% Output:
%   A WIG file containing the genome-wide distribution stored in the
%   cellArray object
%
% Examples:
% write_profile_to_WIG_file(Occ, chrName, 'Occupancy.wig')

noChr = numel(cellArray);
fileID = fopen(WIG_filename, 'w');
for chr = 1:noChr
    fprintf(fileID, 'fixedStep  chrom=%s  start=1  step=1\n', chrName{chr});
    fprintf(fileID, '%0.3f\n', cellArray{chr}(:));
end
fclose(fileID);
