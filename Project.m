%% Figure 3 A (Appendix: 3.2.1 Chemostat growth & 2.5.1 No protein data)
% Project: Improving the phenotype predictions of a yeast genome-scale
% metabolic model by incorporating enzymatic constraints

clear;
clc;
load('ecYeast_v8.6.2.mat')

% find(contains(tmodel.rxnNames,'oxygen')) All the related data
% tmodel.rxnNames(find(contains(tmodel.rxnNames,'oxygen'))) All the name
% find(ismember(tmodel.rxnNames,'oxygen transport via diffusion _extracellular to periplasm')) Find the index of oxygen transport via

id_growth = find(ismember(ecModel.rxnNames,'growth')); % r_2111 
id_biomass = find(ismember(ecModel.rxnNames,'biomass pseudoreaction')); % r_4041 
id_glucose = find(ismember(ecModel.rxnNames,'D-glucose exchange')); % r_1714
id_O2 = find(ismember(ecModel.rxnNames,'oxygen exchange'));  % r_1992
id_CO2 = find(ismember(ecModel.rxnNames,'carbon dioxide exchange')); % r_1672
id_EtOH = find(ismember(ecModel.rxnNames,'ethanol exchange')); % r_1761
id_acetate = find(ismember(ecModel.rxnNames,'acetate exchange')); % r_1634
id_prot_pool_exchange = find(ismember(ecModel.rxns,'prot_pool_exchange'));

%% Under chemostat conditions 5 exchange reactions were limited to physiological levels

id_pyruvate_uptake = find(ismember(ecModel.rxns,'r_2033'));     % 0.05 mmol/gDWh
id_acetate_uptake = find(ismember(ecModel.rxns,'r_1634'));      % 0.62 mmol/gDWh      
id_butanediol_uptake = find(ismember(ecModel.rxns,'r_1549'));       % 1e-5 mmol/gDWh
id_acetaldehyde_uptake = find(ismember(ecModel.rxns,'r_1631'));     % 1e-5 mmol/gDWh
id_glycine_uptake = find(ismember(ecModel.rxns,'r_1810'));          % 1e-5 mmol/gDWh

% ecModel.lb(id_pyruvate_uptake) = 0.05;
% ecModel.lb(id_acetate_uptake) = 0.62;
% ecModel.lb(id_butanediol_uptake) = 1e-5;
% ecModel.lb(id_acetaldehyde_uptake) = 1e-5;
% ecModel.lb(id_glycine_uptake) = 1e-5;

ecModel = setParam(ecModel,'ub','r_2045',0); % Block the L-serine transport between cytoplasm and mitochondria (r_2045)
ecModel = setParam(ecModel,'lb','r_0659',0); % Block the conversion of isocitrate to 2-oxoglutarate in the cytoplasm via NADPH (r_0659)

%% Simulation of a chemostat growth
dilution_rate_range = 0:0.025:0.4;
outV = zeros(numel(ecModel.rxns),numel(dilution_rate_range)); % Creation matrix with zeros
ecModel = setParam(ecModel,'obj','r_1714',1);  % Set the objectove coeff. of glucose 1
totP = -ecModel.lb(id_prot_pool_exchange);

for i = 1:numel(dilution_rate_range)
    model = setParam(ecModel, 'lb','r_2111',dilution_rate_range(i));
    sol=solveLP(model);
    if ~isempty(sol.x)
        model = setParam(model,'lb','r_1714',sol.x(id_glucose)*1.01);
        model = setParam(model,'obj','prot_pool_exchange',1);
        sol=solveLP(model);
        outV(:,i) = sol.x;
    end
end

%% Gather S. cerevisiae experimental data from van Hoeck 1998 
fID = fopen(fullfile(findGECKOroot,'tutorials','full_ecModel','data','vanHoek1998.tsv'),'r');
expData = textscan(fID,'%f %f %f %f %f %f %f %f','Delimiter',';','HeaderLines',2);
fclose(fID);
fluxToPlot = [expData{2} expData{3} expData{4} expData{5} expData{6}]; % O2, CO2, glucose, ethanol, acetate
rxns_indices = [id_O2;id_CO2;id_glucose;id_EtOH; id_acetate];

%% Plot experimental and calculated data
figure
plot(dilution_rate_range,abs(outV(rxns_indices(1),:)),'-',LineWidth=1.5,Color=[0 0.4470 0.7410]);
rectangle('Position', [0.305, 0, max(dilution_rate_range)-0.305, 25], 'FaceColor', [0.85, 0.95, 0.98], 'EdgeColor', 'none');
hold on
plot(dilution_rate_range,abs(outV(rxns_indices(1),:)),'-',LineWidth=1.5,Color=[0 0.4470 0.7410]);
hold on
plot(dilution_rate_range,abs(outV(rxns_indices(2),:)),'-',LineWidth=1.5,Color=[0.4940 0.1840 0.5560]);
hold on
plot(dilution_rate_range,abs(outV(rxns_indices(3),:)),'-',LineWidth=1.5,Color=[0.4660 0.6740 0.1880]);
hold on
plot(dilution_rate_range,abs(outV(rxns_indices(4),:)),'-',LineWidth=1.5,Color=[0.8500 0.3250 0.0980]);
hold on
plot(dilution_rate_range,abs(outV(rxns_indices(5),:)),'-',LineWidth=1.5,Color=[0.9290 0.6940 0.1250]);
hold on
plot(expData{1},fluxToPlot(:,1),markerfacecolor=[0 0.4470 0.7410],Marker='square',LineStyle='none',Color='none')
hold on
plot(expData{1},fluxToPlot(:,2),MarkerFaceColor=[0.4940 0.1840 0.5560],Marker='diamond',LineStyle='none',Color='none')
hold on
plot(expData{1},fluxToPlot(:,3),MarkerFaceColor=[0.4660 0.6740 0.1880],marker='o',LineStyle='none',Color='none')
hold on
plot(expData{1},fluxToPlot(:,4), MarkerFaceColor=[0.8500 0.3250 0.0980],Marker='^',LineStyle='none',Color='none')
hold on
plot(expData{1},fluxToPlot(:,5), markerfacecolor=[0.9290 0.6940 0.1250],Marker='v',LineStyle='none',Color='none')
xline(0.305,'k--',LineWidth=1.5)
legend('O2 consumption','CO2 production','Glucose uptake','Ethanol production','Acetate production','Location','northwest')
legend boxoff
xlabel('Dilution rate [1/h]')
ylim([0 25])
ylabel('Flux [mmol/gDWh]')
hold off
text(0.3, 26.5, sprintf('Only glucose\nlimited'), 'HorizontalAlignment', 'right');
text(0.31, 26.5, sprintf('Glucose and\nprotein limited'));

%% Figure 3 B Predicted top 6 pathways used in terms of mass at increasing dilution rate in aerobic conditions

clear;
clc;
load('ecYeast_v8.6.2.mat')

ecModel = setParam(ecModel, 'obj', 'r_1714', 1);
ecModel = setParam(ecModel, 'lb', params.bioRxn, 0.25);
sol = solveLP(ecModel, 1);
usageData = enzymeUsage(ecModel, sol.x);
usageReport = reportEnzymeUsage(ecModel,usageData);
usageReport.topAbsUsage

