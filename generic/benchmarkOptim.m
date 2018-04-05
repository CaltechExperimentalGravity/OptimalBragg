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
Trials  = 3;

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
    diary('coatingCodeSpeedTestLog.txt')
    settings.Workers          = AllCombos(SimIndex(i),1);
    settings.SwarmSize        = AllCombos(SimIndex(i),2);
    settings.MaxIter          = AllCombos(SimIndex(i),3);
    settings.SelfAdjustment   = AllCombos(SimIndex(i),4);
    settings.SocialAdjustment = AllCombos(SimIndex(i),5);
    settings.TolFun           = AllCombos(SimIndex(i),6);
    
    disp('----------------------------------------------')
    disp(['Iteration ',num2str(i),' of ',num2str(SimNum)])
    disp('Running PSO with the following configuration ...')
    disp(settings)
    disp('----------------------------------------------')
    
    
    for iter = 1:Trials
        OUT = runSwarm_aLIGO_ETM(settings);
        executionTime(iter) = OUT.executionTime;
        fval(iter)          = OUT.fval;
        xout{iter}          = OUT.xout;
    end
    
    optimResults(i).Workers          = settings.Workers;
    optimResults(i).SwarmSize        = settings.SwarmSize;
    optimResults(i).MaxIter          = settings.MaxIter;
    optimResults(i).SelfAdjustment   = settings.SelfAdjustment;
    optimResults(i).SocialAdjustment = settings.SocialAdjustment;
    optimResults(i).TolFun           = settings.TolFun;
    [min_val,min_idx]                = min(fval);
    optimResults(i).minFval          = min_val;
    optimResults(i).minExecutionTime = executionTime(min_idx);
    optimResults(i).OptParams        = xout{min_idx};
    optimResults(i).allExecutionTime = executionTime;
    optimResults(i).allFval          = fval;
    optimResults(i).allParams        = xout;
    
    
    
    % Save Results
    save(sprintf('%s/optimResults.mat',RESULT_FOLDER),'optimResults');
    
    % Rank the results
    T = struct2table(optimResults);
    
    Alpha = 0.7; % Relative Weight (minFval vs minExecutionTime)
    [SCORE,IDX] = sort(Alpha.*normalize(1./(T.minFval),'norm') +  (1-Alpha).*normalize(1./(T.minExecutionTime),'norm'),'descend');
    [~,RANK] = sort(IDX);
    %T2 = table(SCORE,IDX,'VariableNames',{'SCORE','Index'})
    Tfinal = [table(RANK,100*SCORE(RANK),'VariableNames',{'RANK','SCORE'}) T ];
    [~,idx] = sort(Tfinal.RANK,'ascend');
    Tfinal = Tfinal(idx,:);
    writetable(Tfinal,'PSO_benchmark_table.txt');
    disp(Tfinal)
    
    
    % Create HTML table
    Tnew = readtable('PSO_benchmark_table.txt');
    colheads = Tnew.Properties.VariableNames;
    table_cell = [colheads; table2cell(Tnew)];
    caption_str = 'PSO Benchmarking';
    
    html_table(table_cell, 'PSO_benchmark_table.html', 'Caption',caption_str, ...
        'DataFormatStr','%0.2f', 'BackgroundColor','#EFFFFF', 'RowBGColour',{'#000099',[],[],[]}, 'RowFontColour',{'#FFFFB5'}, ...
        'FirstColIsHeading',1);
    
    try
        copyfile('./PSO_benchmark_table.html', '/home/controls/users/public_html/nikhil/');
    catch exception
        disp(getReport(exception))
    end
    
    try
        copyfile('./PSO_benchmark_table.html', RESULT_FOLDER);
    catch exception
        disp(getReport(exception))
    end
    
end

diary('off')
%% (IGNORE)



% % Make Plots (Part of Older code)
% fig111 = figure(111);
% plot(Workers,extractfield(optimResults,'minExecutionTime'),'*-' );
% xlabel('Workers')
% ylabel('Best execution time')
% saveas(gcf,sprintf('./Data/optimResults/%s/optimResults.png',RESULT_FOLDER));
% save(sprintf('./Data/optimResults/%s/optimResults.mat',RESULT_FOLDER),'optimResults');
