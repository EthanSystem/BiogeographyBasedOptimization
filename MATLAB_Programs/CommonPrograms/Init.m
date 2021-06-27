function [InitOPTIONS, MinCost, AvgCost, InitFunction, CostFunction, FeasibleFunction, ...
	MaxParValue, MinParValue, Population, ExportData] = Init(SET, OPTIONS, DisplayFlag, ExportData, ProblemFunction, RandSeed)


% Initialize population-based optimization software.
% WARNING: some of the optimization routines will not work if population size is odd.
% total population size ����������
% generation count limit max number of gene ���ɺ���������� 
% number of genes in each population member   ÿ������Ļ�������
% mutation probability     ͻ�����
%This contains various initialization settings for the optimization methods. 
%You can edit this file to change the population size, the generation count limit, the problem dimension, 
%and the mutation probability of any of the optimization methods that you want to run.


%%
% WARNING: some of the optimization routines will not work if population size is odd.
InitOPTIONS.popsize = OPTIONS.popsize; % total population size ����������
InitOPTIONS.Maxgen = OPTIONS.Maxgen; % generation count limit max number of gene  ���ɺ���������� 
InitOPTIONS.numVar = OPTIONS.numVar; % number of genes in each population member   ÿ������Ļ�������
InitOPTIONS.pmutate = OPTIONS.pmutate; % mutation probability     ͻ�����

%% �����������
if ~exist('RandSeed', 'var')
	RandSeed = round(sum(100*clock));       %���û�ȡʱ��ֵ�����������
end
rand('state', RandSeed); % initialize random number generator  ��ʼ�������������
if DisplayFlag
	disp(['random # seed = ', num2str(RandSeed)]);    
    % ����������Ӵ�����
    ExportData.Details(SET.indexOfOptimizationFunction,SET.indexOfBenchmarkFunction).RandomSeed(SET.monte_index)=RandSeed;


end

%% ��ʼ��
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
MinCost = [Population(1).cost];   % ע�⵽֮ǰ�Ѿ���Population�Ӻõ���������
AvgCost = [AverageCost];

%% ��ʾ�����
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

