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

%params.gR_exp = 0

%% Flux data before NGAM

fluxData = loadFluxData();

bar(fluxData.exchMets, fluxData.exchFluxes)

%% find reactions in models

id_glucose = find(ismember(ecModel.rxnNames,'D-glucose exchange')); % r_1714
id_O2 = find(ismember(ecModel.rxnNames,'oxygen exchange'));  % r_1992
id_CO2 = find(ismember(ecModel.rxnNames,'carbon dioxide exchange')); % r_1672
id_EtOH = find(ismember(ecModel.rxnNames,'ethanol exchange')); % r_1761

id_NGAM = find(ismember(model.rxnNames,'non-growth associated maintenance reaction')); % 'r_4046'
model.ub(id_NGAM) % fixed to 0.7

ecModel.ub(id_O2)

find(contains(model.rxnNames,'associated')); 
model.rxnNames(3414) %NGAM, 'r_4046'

find(contains(model.rxnNames,'biomass'));
a = model.S(3409) %biomass pseudoreaction | 'r_4041'

%% input from GPT how to set GAM

% model.S(:, 3409) = model.S(:, 3409) * (31 / model.c(3409));
% 
% % update growth maximization as the objective function.
% model = setParam(model,'obj','r_4041',1);

%% change GAM function extracted from GECKO
% changeGAM(ecModel,GAM,NGAM)

GAM = 31
bioPos = strcmp(model.rxnNames,'biomass pseudoreaction');
for i = 1:length(model.mets)
    S_ix  = model.S(i,bioPos);
    isGAM = sum(strcmp({'ATP [cytoplasm]','ADP [cytoplasm]','H2O [cytoplasm]', ...
        'H+ [cytoplasm]','phosphate [cytoplasm]'},model.metNames{i})) == 1;
    if S_ix ~= 0 && isGAM
        model.S(i,bioPos) = sign(S_ix)*GAM;
    end
end

%% set the GAM and NGAM
% After the fitting procedure, GAMfitted was equal to 31 mmol/gDW for aerobic conditions
% and 16 mmol/gDW for anaerobic conditions. (from supplementary
% information)

model = setParam(model,'lb','r_4046',19.1);     %max19.1
model = setParam(model,'ub','r_4046',19.1);     %max19.1

ecModel = setParam(ecModel,'lb','r_4046',-0.7);    %max200+
ecModel = setParam(ecModel,'ub','r_4046',0.7);    %max200+
%% experimental data
substrate = ["Glucose";"O2";"CO2";"Ethanol"];
consat38dC = [77.5472161309578;29.4322794340087;-150.905276243908;-122.233091905334];  %mmol/gDW
stdv = [6.1383970094143;6.32511023558826;12.6502204711771;11.1342275574192];
%reference: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4966989/
%from supplementary table 1, %Lahtvee 2016
%table(substrate,consat38dC,stdv) %puts the exp data in a table

%transform consat38dC ??? % not sure how this is done, transforms it about right
%tough for all but the o2 cons.

consat38dC = abs(consat38dC) * 0.1; % mmol/gDW * 1/h
stdv = stdv * 0.1 % mmol/gDW * 1/h

% consat38dC = abs(consat38dC)

%% FBA and manual value readout

% sol = solveLP(ecModel,1);
% bioRxnIdx = getIndexes(ecModel,params.bioRxn,'rxns');
% fprintf('Growth rate: %f /hour.\n', sol.x(bioRxnIdx))

%model
sol = solveLP(model, 1)
printFluxes(model,sol.x, false)
%values @ NGAM 19.1
mod_glucose = -1
mod_O2 = -5.9817
mod_CO2 = 5.9824
mod_EtOH = 0
Yeast = abs([mod_glucose; mod_O2; mod_CO2; mod_EtOH])

%ecModel
sol = solveLP(ecModel)
printFluxes(ecModel,sol.x, false)
%values
ec_glucose = -24.7471/3
ec_O2 = -0.40124
ec_CO2 = 42.9916/3
ec_EtOH = 40.9157/3
ecYeast = abs([ec_glucose; ec_O2; ec_CO2; ec_EtOH])

%% plot

barWidth = 0.3; % Width of each bar

figure;
bar(1:numel(substrate), consat38dC, barWidth, 'y'); % Set barWidth for yellow bars
hold on;


xShift = 0.35; % Shift the bars to the right
xPositionM = (1:numel(substrate)) + xShift;
xPositionE = (1:numel(substrate)) + 2 * xShift;

bar(xPositionM, Yeast, barWidth, 'r'); % Grouped bar for m_values
bar(xPositionE, ecYeast, barWidth, 'b'); % Grouped bar for e_values

errorbar(1:numel(substrate), consat38dC, stdv, '.', 'LineWidth', 0.75, 'Color', 'black');

hold off;

xticks(1:numel(substrate));
xticklabels(substrate);
ylabel('Flux [mmol/gDWh]');
title('growth under temperature stress conditions');

legend('Experimental Data', 'Yeast', 'ecYeast');

