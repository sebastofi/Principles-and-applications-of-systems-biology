%% import models

clear;
clc;

load('Yeast_v8.6.2.mat')
load('ecYeast_v8.6.2.mat')

%% Set ModelAdapter

adapterLocation = fullfile(findGECKOroot,'tutorials','full_ecModel','YeastGEMAdapter.m');
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);

ModelAdapter = ModelAdapterManager.getDefault();
ModelAdapter = params

%% Flux data before NGAM
fluxData = loadFluxData();
bar(fluxData.exchMets, fluxData.exchFluxes) % positive production and negative consumption

%% find reactions in models

id_glucose = find(ismember(ecModel.rxnNames,'D-glucose exchange')); % r_1714
id_O2 = find(ismember(ecModel.rxnNames,'oxygen exchange'));  % r_1992
id_CO2 = find(ismember(ecModel.rxnNames,'carbon dioxide exchange')); % r_1672
id_EtOH = find(ismember(ecModel.rxnNames,'ethanol exchange')); % r_1761

id_NGAM = find(ismember(model.rxnNames,'non-growth associated maintenance reaction')); % 'r_4046'
model.ub(id_NGAM) % fixed to 0.7

ecModel.lb(id_O2)

find(contains(model.rxnNames,'associated')); 
model.rxnNames(3414) %NGAM, 'r_4046'

find(contains(model.rxnNames,'biomass'));
model.S(3409) %biomass pseudoreaction | 'r_4041'

%% set the NGAM
% After the fitting procedure, GAMfitted was equal to 31 mmol/gDW for aerobic conditions
% and 16 mmol/gDW for anaerobic conditions. (from supplementary information)

model = setParam(model,'lb','r_4046',19.1);     %max19.1
model = setParam(model,'ub','r_4046',19.1);     %max19.1

ecModel = setParam(ecModel,'lb','r_4046',0.7);    %max200+  %default=0.7
ecModel = setParam(ecModel,'ub','r_4046',0.7);    %max200+  %default=0.7

%% objective function

ecModel = setParam(ecModel, 'obj', id_glucose, 1);        %set objective function
ecModel = setParam(ecModel, 'lb', params.bioRxn, 0.3);  %growth rate

% model = setParam(model, 'obj', 'r_1714', 1);
% model = setParam(model, 'lb', params.bioRxn, 0.3);

%% experimental data
substrate = ["Glucose";"O2";"CO2";"Ethanol"];
consat38dC = [77.5472161309578;29.4322794340087;-150.905276243908;-122.233091905334];  %mmol/gDW    %what was cons again?
stdv = [6.1383970094143;6.32511023558826;12.6502204711771;11.1342275574192];
%reference: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4966989/
%from supplementary table 1, %Lahtvee 2016
%table(substrate,consat38dC,stdv) %puts the exp data in a table

%transform consat38dC 
consat38dC = abs(consat38dC) * 0.1; % mmol/gDW * 1/h
stdv = stdv * 0.1 % mmol/gDW * 1/h

%% FBA and manual value readout

%model
sol = solveLP(model, 1)
printFluxes(model,sol.x, false)
bioRxnIdx = getIndexes(model,params.bioRxn,'rxns');
fprintf('Growth rate: %f /hour.\n', sol.x(bioRxnIdx))
%values @ NGAM 19.1 | growth rate = 0.000408 (cannot set oj)
mod_glucose = -1
mod_O2 = -5.9817
mod_CO2 = 5.9824
mod_EtOH = 0
Yeast = abs([mod_glucose; mod_O2; mod_CO2; mod_EtOH])
% the Yeast7 model has to have the NGAM at max (19.1) to reach the value on
% the figure, but the growth rate cannot be changed (objective function as
% given by the model itself), so it stays at 0.000408.. this model is built
% probably on top of experimental data, so it has to be at a fixed rate?

%ecModel
sol = solveLP(ecModel)
printFluxes(ecModel,sol.x, false)
bioRxnIdx = getIndexes(ecModel,params.bioRxn,'rxns');
fprintf('Growth rate: %f /hour.\n', sol.x(bioRxnIdx))
%added oj @ NGAM 0.7 | growth rate = 0.3
ec_glucose = -8.3355
ec_O2 = -5.1483
ec_CO2 = 14.3662
ec_EtOH = 6.0755
ecYeast = abs([ec_glucose; ec_O2; ec_CO2; ec_EtOH])
% the ecYeast7 model can have an untouched NGAM of 0.7 in order to reach
% the value on the original figure. However, here we need to set the
% objective function to max, and the growth rate has to be specified to
% 0.1, as the experimental data

%% plot

barWidth = 0.2; % Width of each bar
orange = [1, 0.5, 0];

figure;

xShift = 0.25; % Shift the bars to the right
xPositionM = (1:numel(substrate)) + 2*xShift;
xPositionE = (1:numel(substrate)) + xShift;

bar(1:numel(substrate), consat38dC, barWidth, 'yellow'); % Set barWidth for yellow bars

hold on;
bar(xPositionE, ecYeast, barWidth, 'FaceColor',orange); % Grouped bar for e_values
bar(xPositionM, Yeast, barWidth, 'red'); % Grouped bar for m_values

errorbar(1:numel(substrate), consat38dC, stdv, '.', 'LineWidth', 0.75, 'Color', 'black');
hold off;

xticklabels(substrate);
xticks(1:numel(substrate)); %%%%% label doesnt move this way

ylabel('Flux [mmol/gDWh]');
title('High Energy Demand Fluxes');

legend('Exp. data (0.1 1/h - 38°C)', 'ecYeast8 (0.3 1/h)', 'Yeast8 (0.000408 1/h)');
legend('Location', 'northwest', 'Box', 'off'); % Adjust legend properties

