%% import models

clear;
clc;
load('Yeast_v8.6.2.mat')
load('ecYeast_v8.6.2.mat')
load('ecYeastProt_v8.6.2.mat')

%% Set ModelAdapter

adapterLocation = fullfile(findGECKOroot,'tutorials','full_ecModel','YeastGEMAdapter.m');
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);

ModelAdapter = ModelAdapterManager.getDefault();
%ModelAdapter = params

%% Flux data before NGAM
fluxData = loadFluxData();
%bar(fluxData.exchMets, fluxData.exchFluxes) % positive production and negative consumption

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

%% objective function

ecModel = setParam(ecModel, 'obj', id_glucose, 1);        %set objective function
ecModel = setParam(ecModel, 'lb', params.bioRxn, 0.3);  %growth rate

% As protein model can maximum reach 0.088, also set this as constrain for
% all models.
fluxData.grRate(1) = 0.0880;

% Apply same constraints on exchange fluxes
model = constrainFluxData(model,fluxData,1,'max','loose');
ecModel = constrainFluxData(ecModel,fluxData,1,'max','loose');
ecModelProt = constrainFluxData(ecModelProt,fluxData,1,'max','loose');

solveLP(model)
solveLP(ecModel)
solveLP(ecModelProt)

% Prepare output structure.
minFlux = zeros(numel(model.rxns),3);
maxFlux = minFlux;

% Run ecFVA for each model.
[minFlux(:,1), maxFlux(:,1)] = ecFVA(model, model);
[minFlux(:,2), maxFlux(:,2)] = ecFVA(ecModel, model);
[minFlux(:,3), maxFlux(:,3)] = ecFVA(ecModelProt, model);

% Write results to output file.
output = [model.rxns, model.rxnNames, num2cell([minFlux(:,1), maxFlux(:,1), ...
    minFlux(:,2), maxFlux(:,2), minFlux(:,3), maxFlux(:,3)])]';
fID = fopen(fullfile(params.path,'output','ecFVA.tsv'),'w');
fprintf(fID,'%s %s %s %s %s %s %s %s\n','rxnIDs', 'rxnNames', 'minFlux', ...
            'maxFlux', 'ec-minFlux', 'ec-maxFlux', 'ecP-minFlux', 'ecP-maxFlux');
fprintf(fID,'%s %s %g %g %g %g %g %g\n',output{:});
fclose(fID);


%% Plot ecFVA results and store in output/.
plotEcFVA(minFlux, maxFlux);
saveas(gca, fullfile(params.path,'output','ecFVA.pdf'))


