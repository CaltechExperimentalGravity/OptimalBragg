% Code tp benchmark various optimizers for Coatings Optimization
Workers = [5,10,20,30];
Trials  = 5;

for i = 1:numel(Workers)
    diary('coatingCodeSpeedTestLog.txt')
    settings.Workers = Workers(i);  
    
    for iter = 1:Trials
        OUT = runSwarm_aLIGO_ETM(settings);
        executionTime(iter) = OUT.executionTime;
        fval(iter)          = OUT.fval;
        xout{iter}          = OUT.xout;
    end
    
    optimResults(i).Workers          = Workers(i); 
    [min_val,min_idx]                = min(executionTime);
    optimResults(i).minExecutionTime = min_val;
    optimResults(i).minFval          = fval(min_idx);
    optimResults(i).minParams        = xout{min_idx};
    optimResults(i).allExecutionTime = executionTime;
    optimResults(i).allFval          = fval;
    optimResults(i).allParams        = xout;
    diary('off')
end

% Save Results                                                                                                                                                   
if ~exist('./Data/optimResults')                                                
    mkdir('./Data/optimResults');                                               
end 

DATETIME = char(datetime('now'));
DATETIME = regexprep(DATETIME,':|-| ','_');
RESULT_FOLDER = strcat('optimResults_',DATETIME);
mkdir(sprintf('./Data/optimResults/%s',RESULT_FOLDER));


% Make Plots
fig111 = figure(111);
plot(Workers,extractfield(optimResults,'minExecutionTime'),'*-' );
xlabel('Workers')
ylabel('Best execution time')
saveas(gcf,sprintf('./Data/optimResults/%s/optimResults.png',RESULT_FOLDER));


save(sprintf('./Data/optimResults/%s/optimResults.mat',RESULT_FOLDER),'optimResults');
