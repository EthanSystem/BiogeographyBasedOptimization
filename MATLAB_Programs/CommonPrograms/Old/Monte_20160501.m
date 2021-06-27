function [MeanMinCost, BestMinCost, MeanCPUTime, MeanMinCostNorm, BestMinCostNorm, MeanCPUTimeNorm, ExportData] = Monte(SET , OPTIONS, EXPRESSIONS, ExportData, OptFunction,Bench, DisplayFlag)

% 蒙特卡罗模拟
% Monte Carlo execution of population-based optimization software
% INPUT TimesOfMonte
% 输入值由主函数确定，默认值为函数内部确定的值
% OUTPUT MeanMin is the mean of the best solution found. It is a
% nFunction x nBench array, where nFunction is the number of optimization
% functions that are used, and nBench is the number of benchmarks that
% are optimized.
% OUTPUT MeanMinNorm is MeanMin normalized to a minimum of 1 for each benchmark.
% OUTPUT BestMin is the best solution found by each optimization function
% for each benchmark.
% OUTPUT BestMinNorm is BestMin normalized to a minimum of 1 for each benchmark.
% OUTPUT MeanCPU is the mean CPU time required for each optimization function
% normalized to 1.




% 默认显示输出结果
if ~exist('DisplayFlag','var')
    DisplayFlag = true;
end

nMonte = OPTIONS.TimesOfMonte;    % number of Monte Carlo runs




% compute and display to window
numOfFunction = size(OptFunction,1);
numOfBench = size(Bench,1);
% numOfFunctionOptions = size(OptFunctionOptions, 1);？？？？
MeanMinCost = zeros(numOfFunction, numOfBench);
BestMinCost = inf(numOfFunction, numOfBench);
MeanCPUTime = zeros(numOfFunction, numOfBench);
for indexFunction = 1 : numOfFunction
    for indexBench = 1 : numOfBench
        disp(['Optimization method ', num2str(indexFunction), '/', num2str(numOfFunction), '  , Benchmark function ', num2str(indexBench), '/', num2str(numOfBench)]);
        for indexMonte = 1 : nMonte
            monteTimeStart=tic;
            if DisplayFlag
                disp(['This is the #',num2str(indexMonte),' time']);
            end
            SET.monte_index=indexMonte;
            
            expression01=([OptFunction(indexFunction,:),'(SET,OPTIONS,ExportData,@',Bench(indexBench,:),',',EXPRESSIONS,');']);
            [Cost, ~, ExportData] = eval(expression01);
            
            
            % Calculate cost in each monte
            monteTime=toc(monteTimeStart);
            EachCPUTime(indexMonte)=monteTime;
            EachCost(indexMonte)=Cost(end);
            
        end
        
        % subplot index
        SET.index_of_subplot = SET.index_of_subplot +1;
        % Calculate the mean CPU time , mean min cost , best min cost .
        MeanCPUTime(indexFunction,indexBench)=mean(EachCPUTime);
        MeanMinCost(indexFunction,indexBench)=mean(EachCost);
        BestMinCost(indexFunction,indexBench)=min(EachCost);
        % Export
        ExportData.MeanCPUTime(indexFunction,indexBench)=MeanCPUTime(indexFunction,indexBench);
        ExportData.MeanMinCost(indexFunction,indexBench)=MeanMinCost(indexFunction,indexBench);
        ExportData.BestMinCost(indexFunction,indexBench)=BestMinCost(indexFunction,indexBench);
    
    end
end


% Normalize the results
if min(MeanMinCost) == 0
    MeanMinCostNorm = [];
else
    MeanMinCostNorm = MeanMinCost * diag(1./min(MeanMinCost));
end

if min(BestMinCost) == 0
    BestMinCostNorm = [];
else
    BestMinCostNorm = BestMinCost * diag(1./min(BestMinCost));
end

MeanCPUTime = min(MeanCPUTime);
MeanCPUTimeNorm = MeanCPUTime * diag(1./min(MeanCPUTime));

%Export
ExportData.MeanMinCostNorm = MeanMinCostNorm;
ExportData.BestMinCostNorm = BestMinCostNorm;
ExportData.MeanCPUTimeNorm = MeanCPUTimeNorm;











