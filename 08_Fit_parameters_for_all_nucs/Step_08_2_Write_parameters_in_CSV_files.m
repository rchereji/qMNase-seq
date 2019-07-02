%% S2_exp1
pathName = '../../data/S2_exp1';
load(fullfile(pathName, 'Fit_results.100-200.NF40.S2_exp1.mat'), 'Chr', 'Loc', ...
    'position_std', 'O_vector', 'k1_vector', 'k2_vector');

chrReference = {'2L'; '2R'; '3L'; '3R'; '4'; 'X'; 'Y'; 'rDNA'};
ChrLabel = chrReference(Chr);

StartBP = Loc - 73;
EndBP = Loc + 73; 
T = table(ChrLabel, StartBP, EndBP, position_std, ...
    O_vector, k1_vector, k2_vector);

writetable(T, fullfile(pathName, 'Fit_results.100-200.NF40.S2_exp1.csv'))


%% S2_exp2
pathName = '../../data/S2_exp2';
load(fullfile(pathName, 'Fit_results.100-200.NF35.S2_exp2.mat'), 'Chr', 'Loc', ...
    'position_std', 'O_vector', 'k1_vector', 'k2_vector');

chrReference = {'2L'; '2R'; '3L'; '3R'; '4'; 'X'; 'Y'; 'rDNA'};
ChrLabel = chrReference(Chr);

StartBP = Loc - 73;
EndBP = Loc + 73; 
T = table(ChrLabel, StartBP, EndBP, position_std, ...
    O_vector, k1_vector, k2_vector);

writetable(T, fullfile(pathName, 'Fit_results.100-200.NF35.S2_exp2.csv'))


%% Kc167_exp1
pathName = '../../data/Kc167_exp1';
load(fullfile(pathName, 'Fit_results.100-200.NF35.Kc167_exp1.mat'), 'Chr', 'Loc', ...
    'position_std', 'O_vector', 'k1_vector', 'k2_vector');

chrReference = {'2L'; '2R'; '3L'; '3R'; '4'; 'X'; 'Y'; 'rDNA'};
ChrLabel = chrReference(Chr);

StartBP = Loc - 73;
EndBP = Loc + 73; 
T = table(ChrLabel, StartBP, EndBP, position_std, ...
    O_vector, k1_vector, k2_vector);

writetable(T, fullfile(pathName, 'Fit_results.100-200.NF35.Kc167_exp1.csv'))

