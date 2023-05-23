%%% For random generation of flux seed and plots of optgp sampling (default) with experimental data
model = readCbModel('./iJO1366_paper.xml');

% glucose uptake rate lb setting
model = changeRxnBounds(model, 'EX_glc__D_e', [-15], 'l');

% oxygen uptake flux constraint for aerobic condition
model = changeRxnBounds(model, 'EX_o2_e', [-20], 'l');

%%% to additional biomass flux 
model = changeRxnBounds(model, 'BIOMASS_Ec_iJO1366_WT_53p95M', [0.0], 'b');

model_ori = model;

rng(0, 'twister'); % random number generator initialization

numSample = 1000;


%%%%% Exp data load %%%%%%%%%%%%%%%%%%

% E. coli
Exp_data = dlmread('./Data/Ecoli_exp_data_gur_asr_growth_220222.tsv', '\t', 1, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% random generation of glc uptake, growth, acetate production flux constraints %%%%

C_source_rxn = 'EX_glc__D_e'; % for Ecoli iJO1366

Biomass_rxn = 'BIOMASS_Ec_iJO1366_core_53p95M'; % for Ecoli iJO1366

Target_rxn = 'EX_ac_e'; % for Ecoli iJO1366


% Carbon_source rxn flux random generation
x0 = model.lb(findRxnIDs(model, C_source_rxn)).*rand(numSample, 1);

p = []; % flux pair save matrix

for i = 1:numSample

	
	model = model_ori; % model re-initialization
	
	model = changeRxnBounds(model, C_source_rxn, x0(i), 'b');
	
	fmax = optimizeCbModel(model, 'max');
	
	fmin = optimizeCbModel(model, 'min');
	
	if or(isnan(fmax.f), isnan(fmin.f))
	
		while or(isnan(fmax.f), isnan(fmin.f))
			
			model = model_ori; % model re-initialization
			
			x0(i) = model.lb(findRxnIDs(model, C_source_rxn)).*rand(1);
			
			model = changeRxnBounds(model, C_source_rxn, x0(i), 'b');
			
			fmax = optimizeCbModel(model, 'max');
			
			fmin = optimizeCbModel(model, 'min');
		end
	
	end
	
	
	
	y0 = (fmax.f - fmin.f).*rand(1) + fmin.f;
	
	
	model = changeRxnBounds(model, Biomass_rxn, y0, 'b');
	
	model = changeObjective(model, Target_rxn);
	
	fmax = optimizeCbModel(model, 'max');
	
	fmin = optimizeCbModel(model, 'min');
	
	z0 = (fmax.f - fmin.f).*rand(1) + fmin.f; % effective EX_ac_e flux random generation
	
	temp = [x0(i) y0 z0];
	
	p = [p; temp];
	
end


% data output to csv file
dlmwrite('random_gen_glc_gw_ac_fluxes_for_iJO1366_gur_expansion_ver.csv', p, 'delimiter', ','); % for Ecoli iJO1366 aerobic EX_glc__D_e lb -15 version


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% growth rate vs EX_glc scatter plot

%% With Production Envelope

model = model_ori;

nPts = 100;

x = linspace(model.lb(findRxnIDs(model, C_source_rxn)), 0, nPts);

ymin = zeros(nPts, 1);

ymax = zeros(nPts, 1);


for i = 1:nPts
	
	modelY = changeRxnBounds(model, C_source_rxn, x(i), 'b');
	modelY = changeObjective(modelY, Biomass_rxn);
	
	fmin = optimizeCbModel(modelY, 'min');
	
	fmax = optimizeCbModel(modelY, 'max');
	
	ymin(i) = fmin.f;
	
	ymax(i) = fmax.f;
	
end

% plot random generation flux seeds
scatter(p(:,1), p(:,2));

% plot optgp sample
%scatter(doptgp_sample(:,findRxnIDs(model,C_source_rxn)), doptgp_sample(:,findRxnIDs(model,Biomass_rxn)));

hold on

% plot Exp_data
%s = scatter(Exp_data1(1,:), Exp_data(3,:));
%s.Marker = 'x';
%s.SizeData = 80;

%% Jar
%%%Ishii_et_al
s = scatter(Exp_data(1,1:7), Exp_data(3,1:7))
s.Marker = 'x';
s.SizeData = 80;
%%% Toya
s = scatter(Exp_data(1,10:18), Exp_data(3,10:18))
s.Marker = 'x';
s.SizeData = 80;

%%% Okahashi_et_al
%s = scatter(Exp_data(1,xx:xx), Exp_data(3,xx:xx))
%s.Marker = 'x';
%s.SizeData = 80;


%% Flask
%%% Maeda_et_al
s = scatter(Exp_data(1,8:9), Exp_data(3,8:9))
s.Marker = 'x';
s.SizeData = 80;

%%% Crown_et_al
%s = scatter(Exp_data(1,xx:xx), Exp_data(3,xx:xx))
%s.Marker = 'x';
%s.SizeData = 80;


% plot lines of production envelope
plot(x, ymin, x, ymax, 'LineWidth', 2); % Production Envelope (growth vs glucose uptake) plot


%legend('Default optgp sampling', 'Exp data', 'Minimum Growth', 'Maximum Growth');
legend('Optgp+alpha sampling', 'Jar data1', 'Jar data2', 'Flask data1', 'Minimum Growth', 'Maximum Growth');

xlabel('Specific glucose uptake rate (mmol/gDCW/h)');

ylabel('Specific growth rate (1/h)'); % for Ecoli growth
%ylabel('Specific lactate excretion rate (mmol/gDCW/h)'); % for L. lactis L-lac production

hold off

% for fig file saving as pdf
fig_filename = ['2d_scatter_plot_iJO1366_gw_vs_glc_' num2str(numSample) '_with_production_envelope.pdf']; % for Ecoli iJO1366 aerobic
%fig_filename = ['2d_scatter_plot_iJO1366_anaerobic_ac_vs_gw_' num2str(numSample) '_with_production_envelope.pdf']; % for Ecoli iJO1366 anaerobic
%fig_filename = ['2d_scatter_plot_ecoli_core_gw_vs_glc_default_optgp_' num2str(numSample) '_with_production_envelope_Exp_data.pdf']; % for e_coli_core
%fig_filename = ['2d_scatter_plot_L_lactis_lac_vs_gw_' num2str(numSample) '_with_production_envelope.pdf']; % for L. lactis

saveas(gcf, fig_filename);

% for fig file saving as png
%fig_filename = ['2d_scatter_plot_ecoli_core_gw_vs_glc_default_optgp_' num2str(numSample) '_with_production_envelope_Exp_data.png']; % for e_coli_core
fig_filename = ['2d_scatter_plot_ecoli_iJO1366_gw_vs_glc_optgp_+_alpha' num2str(numSample) '_with_production_envelope_Exp_data.png']; % for E. coli iJO1366

saveas(gcf, fig_filename);

close();






%%%%% EX_ac vs growth rate scatter plot %%%%%%%%%

%% With Production Envelope

model = model_ori;

nPts = 100;

fbasol = optimizeCbModel(model, 'max');

max = fbasol.f;

x = linspace(0, max, nPts);

ymin = zeros(nPts, 1);

ymax = zeros(nPts, 1);


for i = 1:nPts
	
	modelY = changeRxnBounds(model, Biomass_rxn, x(i), 'b');
	modelY = changeObjective(modelY, Target_rxn);
	
	%modelY = changeRxnBounds(model, 'BIOMASS_Ec_iJO1366_core_53p95M', x(i), 'b');
	%modelY = changeObjective(modelY, 'EX_ac_e');
	
	fmin = optimizeCbModel(modelY, 'min');
	
	fmax = optimizeCbModel(modelY, 'max');
	
	ymin(i) = fmin.f;
	
	ymax(i) = fmax.f;
	
end


% plot random generation flux seeds
scatter(p(:,2), p(:,3));

% plot default optgp sample
%scatter(doptgp_sample(:,findRxnIDs(model,Biomass_rxn)), doptgp_sample(:,findRxnIDs(model, Target_rxn)));

hold on

% plot Exp_data
%s = scatter(Exp_data(3,:), Exp_data(2,:));
%s.Marker = 'x';
%s.SizeData = 80;

%% Jar
%%%Ishii_et_al
s = scatter(Exp_data(3,1:7), Exp_data(2,1:7))
s.Marker = 'x';
s.SizeData = 80;

%%% Toya
s = scatter(Exp_data(3,10:18), Exp_data(2,10:18))
s.Marker = 'x';
s.SizeData = 80;

%%% Okahashi_et_al
%s = scatter(Exp_data(3,xx:xx), Exp_data(2,xx:xx))
%s.Marker = 'x';
%s.SizeData = 80;


%% Flask
%%% Maeda_et_al
s = scatter(Exp_data(3,8:9), Exp_data(2,8:9))
s.Marker = 'x';
s.SizeData = 80;

%%% Crown_et_al
%s = scatter(Exp_data(3,xx:xx), Exp_data(2,xx:xx))
%s.Marker = 'x';
%s.SizeData = 80;

% plot lines of production envelope
plot(x, ymin, x, ymax, 'LineWidth', 2); % Production Envelope (Acetate) plot


%legend('Minimum Production', 'Maximum Production', 'Flux Generation');
legend('Optgp+alpha sampling', 'Jar data1', 'Jar data2', 'Flask data1', 'Minimum Production', 'Maximum Production'); % for E. coli iJO1366

xlabel('Specific growth rate (1/h)');

ylabel('Specific acetate excretion rate (mmol/gDCW/h)'); % for Ecoli acetate production
%ylabel('Specific lactate excretion rate (mmol/gDCW/h)'); % for L. lactis L-lac production

hold off

fig_filename = ['2d_scatter_plot_iJO1366_ac_vs_gw_' num2str(numSample) '_with_production_envelope.pdf']; % for Ecoli iJO1366 aerobic
%fig_filename = ['2d_scatter_plot_iJO1366_anaerobic_ac_vs_gw_' num2str(numSample) '_with_production_envelope.pdf']; % for Ecoli iJO1366 anaerobic
%fig_filename = ['2d_scatter_plot_ecoli_core_ac_vs_gw_' num2str(numSample) '_with_production_envelope.pdf']; % for e_coli_core
%fig_filename = ['2d_scatter_plot_L_lactis_lac_vs_gw_' num2str(numSample) '_with_production_envelope.pdf']; % for L. lactis

saveas(gcf, fig_filename);

%% for .png saving
fig_filename = ['2d_scatter_plot_ecoli_iJO1366_ac_vs_gw_optgp+alpha_' num2str(numSample) '_with_production_envelope_Exp_data.png']; % for E. coli iJO1366

saveas(gcf, fig_filename);


close();


