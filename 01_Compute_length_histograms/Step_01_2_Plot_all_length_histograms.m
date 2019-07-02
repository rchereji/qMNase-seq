%% Plot length histograms for S2 cells (rep. 2)

pathName = '../../data/S2_exp2';

f = figure('Position', [50, 50, 450, 250]);
hold all

load(fullfile(pathName, 'Length_histogram.S2_exp2_seq1_1min.dm6.mat'), 'frag_length', 'frag_count')
plot(frag_length, smooth(100*frag_count/sum(frag_count), 5), 'linewidth', 1)

load(fullfile(pathName, 'Length_histogram.S2_exp2_seq1_2min.dm6.mat'), 'frag_length', 'frag_count')
plot(frag_length, smooth(100*frag_count/sum(frag_count), 5), 'linewidth', 1)

load(fullfile(pathName, 'Length_histogram.S2_exp2_seq1_5min.dm6.mat'), 'frag_length', 'frag_count')
plot(frag_length, smooth(100*frag_count/sum(frag_count), 5), 'linewidth', 1)

load(fullfile(pathName, 'Length_histogram.S2_exp2_seq1_15min.dm6.mat'), 'frag_length', 'frag_count')
plot(frag_length, smooth(100*frag_count/sum(frag_count), 5), 'linewidth', 1)

load(fullfile(pathName, 'Length_histogram.S2_exp2_seq1_40min.dm6.mat'), 'frag_length', 'frag_count')
plot(frag_length, smooth(100*frag_count/sum(frag_count), 5), 'linewidth', 1)

load(fullfile(pathName, 'Length_histogram.S2_exp2_seq1_60min.dm6.mat'), 'frag_length', 'frag_count')
plot(frag_length, smooth(100*frag_count/sum(frag_count), 5), 'r', 'linewidth', 1)

xlim([25 250])
ylim([0 1.6])
set(gca, 'FontSize', 11, 'XTick', 0:50:500, 'XGrid', 'on', 'GridLineStyle', '--')
set(gca, 'XMinorTick', 'on', 'TickLength', [0.03 0.025])
xlabel('Fragment length (bp)', 'FontSize', 12)
ylabel('Percentage (%)', 'FontSize', 12)
title('S2 cells')
box off

h = legend({'1 min', '2 min', '5 min', '15 min', '40 min', '60 min'}, 'location', 'EO');
v = get(h,'title');
set(v,'string','Digestion time');
    
set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', 'LengthHistograms_S2_rep2.eps');

%% Plot length histograms for Kc167 cells (rep. 1)

pathName = '../../data/Kc167_exp1';

f = figure('Position', [50, 50, 450, 250]);
hold all

load(fullfile(pathName, 'Length_histogram.Kc167_exp1_seq1_1min.dm6.mat'), 'frag_length', 'frag_count')
plot(frag_length, smooth(100*frag_count/sum(frag_count), 5), 'linewidth', 1)

load(fullfile(pathName, 'Length_histogram.Kc167_exp1_seq1_2min.dm6.mat'), 'frag_length', 'frag_count')
plot(frag_length, smooth(100*frag_count/sum(frag_count), 5), 'linewidth', 1)

load(fullfile(pathName, 'Length_histogram.Kc167_exp1_seq1_5min.dm6.mat'), 'frag_length', 'frag_count')
plot(frag_length, smooth(100*frag_count/sum(frag_count), 5), 'linewidth', 1)

load(fullfile(pathName, 'Length_histogram.Kc167_exp1_seq1_15min.dm6.mat'), 'frag_length', 'frag_count')
plot(frag_length, smooth(100*frag_count/sum(frag_count), 5), 'linewidth', 1)

load(fullfile(pathName, 'Length_histogram.Kc167_exp1_seq1_40min.dm6.mat'), 'frag_length', 'frag_count')
plot(frag_length, smooth(100*frag_count/sum(frag_count), 5), 'linewidth', 1)

load(fullfile(pathName, 'Length_histogram.Kc167_exp1_seq1_60min.dm6.mat'), 'frag_length', 'frag_count')
plot(frag_length, smooth(100*frag_count/sum(frag_count), 5), 'r', 'linewidth', 1)

xlim([25 250])
ylim([0 1.6])
set(gca, 'FontSize', 11, 'XTick', 0:50:500, 'XGrid', 'on', 'GridLineStyle', '--')
set(gca, 'XMinorTick', 'on', 'TickLength', [0.03 0.025])
xlabel('Fragment length (bp)', 'FontSize', 12)
ylabel('Percentage (%)', 'FontSize', 12)
title('Kc167 cells')
box off

h = legend({'1 min', '2 min', '5 min', '15 min', '40 min', '60 min'}, 'location', 'EO');
v = get(h,'title');
set(v,'string','Digestion time');
    
set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', 'LengthHistograms_Kc167_rep2.eps');
