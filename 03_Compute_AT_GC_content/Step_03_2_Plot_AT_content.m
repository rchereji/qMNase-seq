sample_label = {'1min', '2min', '5min', '15min', '40min', '60min'};
noSamples = numel(sample_label);

% Compute the genome-average GC-content
load('AT_GC_content.mat', 'AT_content')
avg_AT_content = nanmean(cat(1, [AT_content{:}]));


%% S2_exp1
figure('Position', [50, 50, 450, 250]);
hold all

pathName = '../../data/S2_exp1';
for s = 1:noSamples
    load(fullfile(pathName, ['AT_GC_distribution.100-200.S2_exp1_', sample_label{s}, '.dm6.mat']), ...
        'ATbinCenters', 'binCounts')
    plot(ATbinCenters, 100*binCounts/sum(binCounts), 'linewidth', 1)
end

set(gca, 'XMinorTick', 'on', 'TickLength', [0.03 0.025])
grid on
set(gca, 'FontSize', 11);
xlim([0.3, 0.85])

xlabel('A/T-content', 'FontSize', 12)
ylabel('Percentage of reads (%)', 'FontSize', 12)
title('S2 cells')

plot([avg_AT_content, avg_AT_content], [0, 10], 'k--')

h = legend({'1 min', '2 min', '5 min', '15 min', '40 min', '60 min'}, 'location', 'EO');
v = get(h,'title');
set(v,'string','Digestion time');
set(gca, 'YGrid', 'on', 'GridLineStyle', '--', 'layer', 'top')

set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(pathName, 'ATcontent_S2_exp1.eps'));


%% S2_exp2
figure('Position', [50, 50, 450, 250]);
hold all

pathName = '../../data/S2_exp2';
for s = 1:noSamples
    load(fullfile(pathName, ['AT_GC_distribution.100-200.S2_exp2_seq1_', sample_label{s}, '.dm6.mat']), ...
        'ATbinCenters', 'binCounts')
    plot(ATbinCenters, 100*binCounts/sum(binCounts), 'linewidth', 1)
end

set(gca, 'XMinorTick', 'on', 'TickLength', [0.03 0.025])
grid on
set(gca, 'FontSize', 11);
xlim([0.3, 0.85])

xlabel('A/T-content', 'FontSize', 12)
ylabel('Percentage of reads (%)', 'FontSize', 12)
title('S2 cells')

plot([avg_AT_content, avg_AT_content], [0, 10], 'k--')

h = legend({'1 min', '2 min', '5 min', '15 min', '40 min', '60 min'}, 'location', 'EO');
v = get(h,'title');
set(v,'string','Digestion time');
set(gca, 'YGrid', 'on', 'GridLineStyle', '--', 'layer', 'top')

set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(pathName, 'ATcontent_S2_exp2.eps'));


%% Kc167_exp1
figure('Position', [50, 50, 450, 250]);
hold all

pathName = '../../data/Kc167_exp1';
for s = 1:noSamples
    load(fullfile(pathName, ['AT_GC_distribution.100-200.Kc167_exp1_seq1_', sample_label{s}, '.dm6.mat']), ...
        'ATbinCenters', 'binCounts')
    plot(ATbinCenters, 100*binCounts/sum(binCounts), 'linewidth', 1)
end

set(gca, 'XMinorTick', 'on', 'TickLength', [0.03 0.025])
grid on
set(gca, 'FontSize', 11);
xlim([0.3, 0.85])

xlabel('A/T-content', 'FontSize', 12)
ylabel('Percentage of reads (%)', 'FontSize', 12)
title('S2 cells')

plot([avg_AT_content, avg_AT_content], [0, 10], 'k--')

h = legend({'1 min', '2 min', '5 min', '15 min', '40 min', '60 min'}, 'location', 'EO');
v = get(h,'title');
set(v,'string','Digestion time');
set(gca, 'YGrid', 'on', 'GridLineStyle', '--', 'layer', 'top')

set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', fullfile(pathName, 'ATcontent_Kc167_exp1.eps'));
