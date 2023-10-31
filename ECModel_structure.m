%% STAGE 1: expansion from a starting metabolic model to an ecModel structure
clc
clear

% Here we are making a new ecModel structure becuase the protocol.m given
% on github is performed with the new version of GECKO (v. 3), whereas the
% ecModel strcuture from the research article was processed with GECKO (v.
% 1). Therefore, the ecModel from the article is not compatible with GECKO
% 3. 

% STEP 1 Set modelAdapter
adapterLocation = fullfile(findGECKOroot,'tutorials','full_ecModel','YeastGEMAdapter.m');
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);

ModelAdapter = ModelAdapterManager.getDefault();
params = ModelAdapter.getParameters();

% STEP 2 Load conventional yeast-GEM
model = loadConventionalGEM();

% STEP 3 Prepare ecModel
[ecModel, noUniprot] = makeEcModel(model,false);

% STEP 4 Annotate with complex data
[ecModel, foundComplex, proposedComplex] = applyComplexData(ecModel);

%% STAGE 2: integration of kcat into the ecModel structure
% STEP 6 Fuzzy matching with BRENDA
ecModel         = getECfromGEM(ecModel);
noEC = cellfun(@isempty, ecModel.ec.eccodes);
ecModel         = getECfromDatabase(ecModel,noEC);
ecModel         = getECfromDatabase(ecModel);
kcatList_fuzzy  = fuzzyKcatMatching(ecModel);

% STEP 7 DLKcat prediction through machine learning
[ecModel, noSmiles] = findMetSmiles(ecModel);
kcatList_DLKcat = readDLKcatOutput(ecModel);

% STEP 8 Combine kcat from BRENDA and DLKcat
kcatList_merged = mergeDLKcatAndFuzzyKcats(kcatList_DLKcat, kcatList_fuzzy);

% STEP 9 Take kcatList and populate edModel.ec.kcat
ecModel  = selectKcatValue(ecModel, kcatList_merged);

% STEP 10 Apply custom kcat values
[ecModel, rxnUpdated, notMatch] = applyCustomKcats(ecModel);

% STEP 11 Get kcat values across isozymes
ecModel = getKcatAcrossIsozymes(ecModel);

% STEP 12 Get standard kcat
[ecModel, rxnsMissingGPR, standardMW, standardKcat] = getStandardKcat(ecModel);

% SAVE
filename = 'ecYeast_v8.6.2.mat';
save(filename, 'ecModel','params');

%%   This part was just for testing
ecModel = applyKcatConstraints(ecModel);
Ptot  = params.Ptot;
f     = params.f;
sigma = params.sigma;
ecModel = setProtPoolSize(ecModel,Ptot,f,sigma);


id_dilution_rate = find(ismember(ecModel.rxnNames,'biomass pseudoreaction'));      
id_glucose_uptake = find(ismember(ecModel.rxnNames,'D-glucose exchange'));
ecModel = setParam(ecModel,'lb','r_1714',-1000);
ecModel = setParam(ecModel,'obj','r_4041',1);
sol = solveLP(ecModel,1);
bioRxnIdx = getIndexes(ecModel,params.bioRxn,'rxns');
fprintf('Growth rate: %f /hour.\n', sol.x(bioRxnIdx))


