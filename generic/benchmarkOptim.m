% Code to find the best hyperparameters for PSO based Coatings Optimization
% The following code grid the n-dimensional parameter space & 
% perform brute force search over randomly selected indices
% Nikhil (5th April 2018)

% PSO Parameters
Workers = [5,20,40];
Workers = [0 Workers];
SwarmSize = [5,50,500]; % times nvars
MaxIter = [100,500,1000];
SelfAdjustment = [0.1,1,10];
SocialAdjustment = [0.1,1,10];
TolFun  = [1e-8,1e-6,1e-3];

% Generate All Possible Combos for Grid Search
AllCombos = allcomb(Workers,SwarmSize,MaxIter,SelfAdjustment,...
          SocialAdjustment,TolFun);

% Number of trials       
Trials  = 1;

% Generate random indices to search
SimNum = 100; % Total number ofd simulations
rng(0,'twister');  
SimIndex =  randi(size(AllCombos,1),SimNum,1);

% Create Result Folder                                                                                                                                                   
if ~exist('./Data/benchmark_Results','dir')                                                
    mkdir('./Data/benchmark_Results');                                               
end 

DATETIME = char(datetime('now'));
DATETIME = regexprep(DATETIME,':|-| ','_');
RESULT_FOLDER = strcat('benchmark_Results_',DATETIME);
RESULT_FOLDER = sprintf('./Data/benchmark_Results/%s',RESULT_FOLDER);
mkdir(RESULT_FOLDER);

%% Start the Grid Run
for i = 1:SimNum
    %diary('coatingCodeSpeedTestLog.txt')
    settings.Workers          = AllCombos(SimIndex(i),1);
    settings.SwarmSize        = AllCombos(SimIndex(i),2);
    settings.MaxIter          = AllCombos(SimIndex(i),3);
    settings.SelfAdjustment   = AllCombos(SimIndex(i),4);
    settings.SocialAdjustment = AllCombos(SimIndex(i),5);
    settings.TolFun           = AllCombos(SimIndex(i),6);
    
    
    for iter = 1:Trials
        OUT = runSwarm_aLIGO_ETM(settings);
        executionTime(iter) = OUT.executionTime;
        fval(iter)          = OUT.fval;
        xout{iter}          = OUT.xout;
    end
    
    optimResults(i).Workers          = settings.Workers; 
    [min_val,min_idx]                = min(fval);
    optimResults(i).minFval          = min_val;
    optimResults(i).minExecutionTime = executionTime(min_idx);
    optimResults(i).OptParams        = xout{min_idx};
    optimResults(i).allExecutionTime = executionTime;
    optimResults(i).allFval          = fval;
    optimResults(i).allParams        = xout;
    %diary('off')
    
    % Save Results
    save(sprintf('%s/optimResults.mat',RESULT_FOLDER),'optimResults');
    
end




% % Make Plots (Part of Older code)
% fig111 = figure(111);
% plot(Workers,extractfield(optimResults,'minExecutionTime'),'*-' );
% xlabel('Workers')
% ylabel('Best execution time')
% saveas(gcf,sprintf('./Data/optimResults/%s/optimResults.png',RESULT_FOLDER));
% save(sprintf('./Data/optimResults/%s/optimResults.mat',RESULT_FOLDER),'optimResults');
