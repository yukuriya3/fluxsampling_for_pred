%%% Getting Correlated Reaction Sets

initCobraToolbox(false);

model = readCbModel('../e_core_core.xml');

samples = load('Ecoli_core_merged_sample.csv');


%[setsSorted, setsNoSorted, setSize] = identifyCorrelSets(model, samples', 0.9025);
[setsSorted, setsNoSorted, setSize] = identifyCorrelSets(model, samples', 0.810);

setNames = [];

setNumbers = [];

disp('Correlated Reaction Sets');

for i = 1:length(setsSorted)
	setNames = [setNames; setsSorted{i}.names];
	setNumbers = [setNumbers; i * ones(length(setsSorted{i}.names), 1)];
	
	disp([i, setsSorted{i}.names']);
end



