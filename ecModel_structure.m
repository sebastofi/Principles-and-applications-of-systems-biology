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

% STEP 13 Apply kcat constraints from ecModel.ec.kcat to ecModel.S
ecModel = applyKcatConstraints(ecModel);

% STEP 14 Set upper bound of protein pool
params.Ptot =  0.4005; % g/gDW
params.f = 0.4461; % g/g
params.sigma = 0.51; % 51 per-cent saturation
ecModel = setProtPoolSize(ecModel,params.Ptot,params.f,params.sigma);

%% STAGE 3: model tuning
% STEP 15 Test maximum growth rate
ecModel = setParam(ecModel,'lb','r_1714',-1000);
ecModel = setParam(ecModel,'obj','r_4041',1);
sol = solveLP(ecModel,1);
bioRxnIdx = getIndexes(ecModel,params.bioRxn,'rxns');
fprintf('Growth rate: %f /hour.\n', sol.x(bioRxnIdx))

% STEP 17 Sensitivity tuning
ecModel = setProtPoolSize(ecModel);
[ecModel, tunedKcats] = sensitivityTuning(ecModel);
struct2table(tunedKcats)

% STEP 18 Curate kcat values based on kcat tuning
rxnIdx = find(strcmp(kcatList_merged.rxns,'r_0079'));
doc fuzzyKcatMatching % To check the meaning of wildcardLvl and origin.
kcatList_merged.wildcardLvl(rxnIdx) % 0: no EC number wildcard.
kcatList_merged.origin(rxnIdx) % 4: any organism, any substrate, kcat.
kcatList_merged.eccodes(rxnIdx) % EC number 6.3.5.3.

enzMW = ecModel.ec.mw(strcmp(ecModel.ec.enzymes,'P38972')); % Get MW of the enzyme.
convKcat = 2.15; % umol/min/mg protein, same as mmol/min/g protein.
convKcat = convKcat / 1000; % mol/min/g protein.
convKcat = convKcat / 60; % mol/sec/g protein.
convKcat = convKcat * enzMW; % mol/sec/mol protein, same as 1/sec.
ecModel = setKcatForReactions(ecModel,'r_0079',convKcat);
ecModel = applyKcatConstraints(ecModel);

% SAVE
filename = 'Yeast_v8.6.2.mat';
save (filename,'model')
filename = 'ecYeast_v8.6.2.mat';
save(filename, 'ecModel','params');