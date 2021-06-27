function [InitOPTIONS, MinCost, AvgCost, InitFunction, CostFunction, FeasibleFunction, ...
	MaxParValue, MinParValue, Population, ExportData] = Init(SET, OPTIONS, DisplayFlag, ExportData, ProblemFunction, RandSeed)


% Initialize population-based optimization software.
% WARNING: some of the optimization routines will not work if population size is odd.
% total population size 个体总数量
% generation count limit max number of gene 生成后代代数限制 
% number of genes in each population member   每个个体的基因数量
% mutation probability     突变概率
%This contains various initialization settings for the optimization methods. 
%You can edit this file to change the population size, the generation count limit, the problem dimension, 
%and the mutation probability of any of the optimization methods that you want to run.


%%
% WARNING: some of the optimization routines will not work if population size is odd.
InitOPTIONS.popsize = OPTIONS.popsize; % total population size 个体总数量
InitOPTIONS.Maxgen = OPTIONS.Maxgen; % generation count limit max number of gene  生成后代代数限制 
InitOPTIONS.numVar = OPTIONS.numVar; % number of genes in each population member   每个个体的基因数量
InitOPTIONS.pmutate = OPTIONS.pmutate; % mutation probability     突变概率

%% 产生随机种子
if ~exist('RandSeed', 'var')
	RandSeed = round(sum(100*clock));       %利用获取时间值生成随机种子
end
rand('state', RandSeed); % initialize random number generator  初始化随机数生成器
if DisplayFlag
	disp(['random # seed = ', num2str(RandSeed)]);    
    % 把随机数种子存起来
    ExportData.Details(SET.indexOfOptimizationFunction,SET.indexOfBenchmarkFunction).RandomSeed(SET.monte_index)=RandSeed;


end

%% 初始化
% Get the addresses of the initialization, cost, and feasibility functions.
[InitFunction, CostFunction, FeasibleFunction] = ProblemFunction();
% Initialize the population.
[MaxParValue, MinParValue, Population, InitOPTIONS] = InitFunction(InitOPTIONS);
% Make sure the population does not have duplicates.
Population = ClearDups(Population, MaxParValue, MinParValue);
% Compute cost of each individual
Population = CostFunction(InitOPTIONS, Population);
% Sort the population from most fit to least fit
Population = PopSort(Population);
% Compute the average cost
AverageCost = ComputeAveCost(Population);
% Display info to screen
MinCost = [Population(1).cost];   % 注意到之前已经将Population从好到坏排序了
AvgCost = [AverageCost];

%% 显示与输出
if DisplayFlag
	disp(['The best and mean of Generation # 0 are ', num2str(MinCost(end)), ' and ', num2str(AvgCost(end))]);
end

for i=1:InitOPTIONS.popsize
    ExportData.Details(SET.indexOfOptimizationFunction,SET.indexOfBenchmarkFunction).EachPopulationInEachGenerationAndMonte(1,SET.monte_index).Population(i).chrom=Population(i).chrom;
    ExportData.Details(SET.indexOfOptimizationFunction,SET.indexOfBenchmarkFunction).EachPopulationInEachGenerationAndMonte(1,SET.monte_index).Population(i).FitnessCost=0;
    ExportData.Details(SET.indexOfOptimizationFunction,SET.indexOfBenchmarkFunction).EachPopulationInEachGenerationAndMonte(1,SET.monte_index).Population(i).Cost=Population(i).cost;
end
ExportData.Details(SET.indexOfOptimizationFunction,SET.indexOfBenchmarkFunction).BestAndMeanInEachGenerationAndMonte(1).generation(SET.monte_index)=0;
ExportData.Details(SET.indexOfOptimizationFunction,SET.indexOfBenchmarkFunction).BestAndMeanInEachGenerationAndMonte(1).best(SET.monte_index)=MinCost(end);
ExportData.Details(SET.indexOfOptimizationFunction,SET.indexOfBenchmarkFunction).BestAndMeanInEachGenerationAndMonte(1).mean(SET.monte_index)=AvgCost(end);




return;

