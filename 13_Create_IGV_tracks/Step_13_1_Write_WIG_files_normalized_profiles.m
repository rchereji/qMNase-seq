% S2_exp1
pathName = '../../data/S2_exp1';
allFiles = dir([pathName, '/Spikein_normalized_profiles.*.mat']);

for f = 1:length(allFiles)
    filename = fullfile(pathName, allFiles(f).name);
    load(filename, 'Dyads', 'Occ_101', 'chrName')
    
    [~, ind, ~] = intersect(chrName, {'2L', '2R', '3L', '3R', '4', 'X', 'Y'}, 'stable');
    Dyads = Dyads(ind);
    Occ_101 = Occ_101(ind);
    chrName = chrName(ind);
    chrName = cellfun(@(s) ['chr', s], chrName, 'UniformOutput', 0);
    
    setName = allFiles(f).name(29:end-4);
    write_profile_to_WIG_file(Dyads, chrName, fullfile(pathName, ['Dyads.', setName, '.wig']));
    write_profile_to_WIG_file(Occ_101, chrName, fullfile(pathName, ['Occ_101bp.', setName, '.wig']));
    disp(f)
end


%% S2_exp2
pathName = '../../data/S2_exp2';
allFiles = dir([pathName, '/Spikein_normalized_profiles.*.mat']);

for f = 1:length(allFiles)
    filename = fullfile(pathName, allFiles(f).name);
    load(filename, 'Dyads', 'Occ_101', 'chrName')
    
    [~, ind, ~] = intersect(chrName, {'2L', '2R', '3L', '3R', '4', 'X', 'Y'}, 'stable');
    Dyads = Dyads(ind);
    Occ_101 = Occ_101(ind);
    chrName = chrName(ind);
    chrName = cellfun(@(s) ['chr', s], chrName, 'UniformOutput', 0);
    
    setName = allFiles(f).name(29:end-4);
    write_profile_to_WIG_file(Dyads, chrName, fullfile(pathName, ['Dyads.', setName, '.wig']));
    write_profile_to_WIG_file(Occ_101, chrName, fullfile(pathName, ['Occ_101bp.', setName, '.wig']));
    disp(f)
end


%% Kc167_exp1
pathName = '../../data/Kc167_exp1';
allFiles = dir([pathName, '/Spikein_normalized_profiles.*.mat']);

for f = 1:length(allFiles)
    filename = fullfile(pathName, allFiles(f).name);
    load(filename, 'Dyads', 'Occ_101', 'chrName')
    
    [~, ind, ~] = intersect(chrName, {'2L', '2R', '3L', '3R', '4', 'X', 'Y'}, 'stable');
    Dyads = Dyads(ind);
    Occ_101 = Occ_101(ind);
    chrName = chrName(ind);
    chrName = cellfun(@(s) ['chr', s], chrName, 'UniformOutput', 0);
    
    setName = allFiles(f).name(29:end-4);
    write_profile_to_WIG_file(Dyads, chrName, fullfile(pathName, ['Dyads.', setName, '.wig']));
    write_profile_to_WIG_file(Occ_101, chrName, fullfile(pathName, ['Occ_101bp.', setName, '.wig']));
    disp(f)
end
