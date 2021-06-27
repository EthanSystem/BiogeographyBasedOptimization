function [MinCost, Hamming] = BBO(ProblemFunction, DisplayFlag, ProbFlag, RandSeed)

% Biogeography-based optimization (BBO) software for minimizing a general function
% BBO算法用于最小化通用函数的解的误差值

% INPUTS: ProblemFunction is the handle of the function that returns
%         the handles of the initialization, cost, and feasibility
%         functions.  初始化、价值、可行性函数的句柄
%         DisplayFlag = true or false, whether or not to display and plot results.
%         ProbFlag = true or false, whether or not to use probabilities to update emigration rates.
%         RandSeed = random number seed
% OUTPUTS: MinCost = array of best solution, one element for each generation
%          Hamming = final Hamming distance between solutions
%          解与解之间的Hamming距离
% CAVEAT: The "ClearDups" function that is called below replaces duplicates with randomly-generated
%         individuals, but it does not then recalculate the cost of the replaced individuals.


if ~exist('DisplayFlagpar', 'var')
	DisplayFlag = true;
end
if ~exist('ProbFlag', 'var')
	ProbFlag = false;
end
if ~exist('RandSeed', 'var')
	RandSeed = round(sum(100*clock));
end

[OPTIONS, MinCost, AvgCost, InitFunction, CostFunction, FeasibleFunction, ...
	MaxParValue, MinParValue, Population] = Init(DisplayFlag, ProblemFunction, RandSeed);

Population = CostFunction(OPTIONS, Population);

OPTIONS.pmodify = 1; % habitat modification probability 栖息地的变动概率

OPTIONS.pmutate = 0.005; % initial mutation probability初始化突变概率

Keep = 2; % elitism parameter: how many of the best habitats to keep from one generation to the next 精英参数
lambdaLower = 0.0; % lower bound for immigration probabilty per gene
lambdaUpper = 1.0; % upper bound for immigration probabilty per gene
dt = 1; % step size used for numerical integration of probabilities
I = 1; % max immigration rate for each island
E = 1; % max emigration rate, for each island
P = OPTIONS.popsize; % max species count, for each island


% Initialize the species count probability of each habitat
% Later we might want to initialize probabilities based on cost
%   初始化每一个栖息地的种群数量的概率
for j = 1 : length(Population)
	Prob(j) = 1 / length(Population);
end



        
%% Begin the optimization loop
for GenIndex = 1 : OPTIONS.Maxgen   % 从第一代开始
	%% Save the best habitats in a temporary array.
	for j = 1 : Keep
		chromKeep(j,:) = Population(j).chrom;
		costKeep(j) = Population(j).cost;
    end
    
	%% Map cost values to species counts.
	[Population] = GetSpeciesCounts(Population, P);
    
	%% Compute immigration rate and emigration rate for each species count.
	% lambda(i) is the immigration rate for habitat i.
	% mu(i) is the emigration rate for habitat i.
	[lambda, mu] = GetLambdaMu(Population, I, E, P);
    
     %%	 Compute the time derivative of Prob(i) for each habitat i.
		% 	计算每个栖息地的概率对时间的导数
        % 注意：再次强调，运行函数前，已经假定栖息地 i已经从优到劣排序好了。
	if ProbFlag         % If allowed to use probablities to update emigration rates . 
		for j = 1 : length(Population)
			% Compute lambda for one less than the species count of habitat i.
			lambdaMinus = I * (1 - (Population(j).SpeciesCount - 1) / P);   
			% Compute mu for one more than the species count of habitat i.
			muPlus = E * (Population(j).SpeciesCount + 1) / P;
			% Compute Prob for one less than and one more than the species count of habitat i.
			% Note that species counts are arranged in an order opposite to that presented in
			% MacArthur and Wilson's book - that is, the most fit
			% habitat has index 1, which has the highest species count.
			if j < length(Population)
				ProbMinus = Prob(j+1);
			else
				ProbMinus = 0;
			end
			if j > 1
				ProbPlus = Prob(j-1);
            else
				ProbPlus = 0;
			end
            
			ProbDot(j) = -(lambda(j) + mu(j)) * Prob(j) + lambdaMinus * ProbMinus + muPlus * ProbPlus;
        end
        
		% Compute the new probabilities for each species count.
		Prob = Prob + ProbDot * dt;
		Prob = max(Prob, 0);    % 使Prob的所有值非负
		Prob = Prob / sum(Prob);                
    end
    
    
    
	%% Now use lambda and mu to decide how much information to share between habitats.
    % 迁移运算：迁入主导部分变量迁移。用 lambda 和 mu 确定在 habitats 之间可以分享的信息量。
	lambdaMin = min(lambda);
	lambdaMax = max(lambda);
	for k = 1 : length(Population)
		if rand > OPTIONS.pmodify
			continue;
		end
		% Normalize the immigration rate.      
        % lamdaScale 将现有的 lambda(k) 的尺度缩放到 [0,1] 实数空间，使之可以与rand()比较
		lambdaScale = lambdaLower + (lambdaUpper - lambdaLower) * (lambda(k) - lambdaMin) / (lambdaMax - lambdaMin);
		%% Probabilistically input new information into habitat i
		for j = 1 : OPTIONS.numVar
			if rand < lambdaScale    % 根据lambdaScale确定这个特征值是否迁徙
				% Pick a habitat from which to obtain a feature
                %% 根据mu 选取一个 habitat ，通过轮盘赌选取法，mu越大的habitat越容易被选中
				RandomNum = rand * sum(mu);
				Select = mu(1);
				SelectIndex = 1;
				while (Select < RandomNum) & (SelectIndex < OPTIONS.popsize)
					SelectIndex = SelectIndex + 1;
					Select = Select + mu(SelectIndex);
                end
               %发生迁徙，接下来将把选中的栖息地的第 j 个SIV复制给第k个栖息地
				Island(k,j) = Population(SelectIndex).chrom(j);     
			else
				Island(k,j) = Population(k).chrom(j);       % 不发生迁徙
			end
		end
    end
    
    %%  Mutation 设置 突变
	if ProbFlag
		% Mutation
		Pmax = max(Prob);
		MutationRate = OPTIONS.pmutate * (1 - Prob / Pmax);
		% Mutate only the worst half of the solutions
		Population = PopSort(Population);
		for k = round(length(Population)/2) : length(Population)    % 此时Population已经排好序了
			for parnum = 1 : OPTIONS.numVar
				if rand < MutationRate(k)
					Island(k,parnum) = floor(MinParValue + (MaxParValue - MinParValue + 1) * rand);  
				end
			end
		end
    end
    
    %% 更新
	% Replace the habitats with their new versions.
    for k = 1 : length(Population)
		Population(k).chrom = Island(k,:);
    end
	% Make sure each individual is legal.
	Population = FeasibleFunction(OPTIONS, Population);
	% Calculate cost
	Population = CostFunction(OPTIONS, Population);
	% Sort from best to worst
	Population = PopSort(Population);
	% Replace the worst with the previous generation's elites.
	n = length(Population);
	for k = 1 : Keep
		Population(n-k+1).chrom = chromKeep(k,:);
		Population(n-k+1).cost = costKeep(k);
	end
	% Make sure the population does not have duplicates.
	Population = ClearDups(Population, MaxParValue, MinParValue);
	% Sort from best to worst
	Population = PopSort(Population);
	% Compute the average cost
	[AverageCost, nLegal] = ComputeAveCost(Population);
	% Display info to screen
	MinCost = [MinCost Population(1).cost];
	AvgCost = [AvgCost AverageCost];
	if DisplayFlag
		disp(['The best and mean of Generation # ', num2str(GenIndex), ' are ',...
			num2str(MinCost(end)), ' and ', num2str(AvgCost(end))]);
	end
end

%% 结论
Conclude(DisplayFlag, OPTIONS, Population, nLegal, MinCost);
% Obtain a measure of population diversity
% 从栖息地k到下一层级基因的遍历
for k = 1 : length(Population)
	Chrom = Population(k).chrom;
	for j = MinParValue : MaxParValue
		indices = find(Chrom == j);
		CountArr(k,j) = length(indices); % array containing gene counts of each habitat
	end
end
Hamming = 0;
for m = 1 : length(Population)
	for j = m+1 : length(Population)
		for k = MinParValue : MaxParValue
			Hamming = Hamming + abs(CountArr(m,k) - CountArr(j,k));
		end
	end
end
if DisplayFlag
	disp(['Diversity measure = ', num2str(Hamming)]);
end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%GetSpecies函数
function [Population] = GetSpeciesCounts(Population, P)

% Map cost values to species counts.
% This loop assumes the population is already sorted from most fit to least fit.
for i = 1 : length(Population)
	if Population(i).cost < inf
		Population(i).SpeciesCount = P - i;
	else
		Population(i).SpeciesCount = 0;
	end
end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%GetLambdaMu函数
function [lambda, mu] = GetLambdaMu(Population, I, E, P)

% Compute immigration rate and extinction rate for each species count.
% lambda(i) is the immigration rate for individual i.
% mu(i) is the extinction rate for individual i.

for i = 1 : length(Population)
	lambda(i) = I * (1 - Population(i).SpeciesCount / P);
	mu(i) = E * Population(i).SpeciesCount / P;
end
return;
