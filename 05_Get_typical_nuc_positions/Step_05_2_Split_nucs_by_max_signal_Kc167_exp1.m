% Load nuc. counts
pathName = '../../data/Kc167_exp1';
load(fullfile(pathName, 'NucCounts.100-200.Kc167_exp1.mat'), 'Chr', 'Loc', 'NucCount')

% Get the timepoint of max. nuc. count for each nucleosome
[~, max_ids] = arrayfun(@(i) max(NucCount(i,:)), 1:size(NucCount, 1));


%% Plot the distribution of nuc. types
figure('Position', [50, 50, 400, 300])
edges = 0.5:1:6.5;
Counts = histc(max_ids, edges); Counts = Counts(1:6);
y = 100*Counts/sum(Counts);
bar(1:6, y, 1, 'FaceColor', [.85 .85 .85]);

xlim([0.5, 6.5])

set(gca, 'FontSize', 12, 'xticklabel', {'1', '2', '5', '15', '40', '60'})
xlabel('t_{max} (minutes)', 'FontSize', 14)
ylabel('Percentage of nucleosomes (%)', 'FontSize', 14)
title('Kc167_exp1', 'interpreter', 'none')

set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', ...
    fullfile(pathName, 'Percentage_of_color-coded_nucs.Kc167_exp1.eps'));


%% Save results to a BED file
chrReference = {'2L', '2R', '3L', '3R', '4', 'X', 'Y', 'rDNA'};
nucColors = {'230,25,75'; ...
             '245,130,48'; ...
             '255,225,25'; ...
             '210,245,60'; ...
             '60,180,75'; ...
             '0,130,200'};
chrLabel = chrReference(Chr)';

noNucs = numel(Chr);

fileID = fopen(fullfile(pathName, 'Color_coded_nucs.Kc167_exp1.bed'), 'W');
header_txt = 'track name="Colored nucs" description="Color-coded nucleosome positions" useScore=0 itemRgb="On"';
fprintf(fileID, '%s\n', header_txt);
Fmt = '%s\t%d\t%d\tNuc_%d\t0\t*\t%d\t%d\t%s\n';
for n = 1:noNucs
  fprintf(fileID, Fmt, ...
      chrLabel{n}, Loc(n)-74, Loc(n)+73, n, Loc(n)-74, Loc(n)+73, nucColors{max_ids(n)});
end
fclose(fileID);
