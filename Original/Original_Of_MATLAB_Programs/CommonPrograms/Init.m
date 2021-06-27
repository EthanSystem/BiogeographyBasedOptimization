function [OPTIONS, MinCost, AvgCost, InitFunction, CostFunction, FeasibleFunction, ...
	MaxParValue, MinParValue, Population] = Init(DisplayFlag, ProblemFunction, RandSeed)

% Initialize population-based optimization software.
% WARNING: some of the optimization routines will not work if population size is odd.
% total population size ����������
% generation count limit max number of gene ���ɺ���������� 
% number of genes in each population member   ÿ������Ļ�������
% mutation probability     ͻ�����
%This contains various initialization settings for the optimization methods. 
%You can edit this file to change the population size, the generation count limit, the problem dimension, 
%and the mutation probability of any of the optimization methods that you want to run.


% WARNING: some of the optimization routines will not work if population size is odd.
OPTIONS.popsize = 10; % total population size ����������
OPTIONS.Maxgen = 1000; % generation count limit   ���ɺ���������� max number of gene
OPTIONS.numVar = 3; % number of genes in each population member   ÿ������Ļ�������
OPTIONS.pmutate = 0.05; % mutation probability     ͻ�����

if ~exist('RandSeed', 'var')
	RandSeed = round(sum(100*clock));       %���û�ȡʱ��ֵ�����������
end
rand('state', RandSeed); % initialize random number generator  ��ʼ�������������
if DisplayFlag
	disp(['random # seed = ', num2str(RandSeed)]);    
end

% Get the addresses of the initialization, cost, and feasibility functions.
% ����ʽ���ұ߸ĳ���Ҫ���Եĺ���
[InitFunction, CostFunction, FeasibleFunction] = ProblemFunction();
% Initialize the population.
[MaxParValue, MinParValue, Population, OPTIONS] = InitFunction(OPTIONS);
% Make sure the population does not have duplicates.
Population = ClearDups(Population, MaxParValue, MinParValue);
% Compute cost of each individual
Population = CostFunction(OPTIONS, Population);
% Sort the population from most fit to least fit
Population = PopSort(Population);
% Compute the average cost
AverageCost = ComputeAveCost(Population);
% Display info to screen
MinCost = [Population(1).cost];
AvgCost = [AverageCost];
if DisplayFlag
	disp(['The best and mean of Generation # 0 are ', num2str(MinCost(end)), ' and ', num2str(AvgCost(end))]);
end

return;
