function [MinCost] = ES(ProblemFunction, DisplayFlag)

% (mu-plus-lambda) Evolutionary Strategy for optimizing a general function.

% INPUTS: ProblemFunction is the handle of the function that returns 
%         the handles of the initialization, cost, and feasibility functions.
%         DisplayFlag says whether or not to display information during iterations and plot results.

if ~exist('DisplayFlag', 'var')
    DisplayFlag = true;
end

[OPTIONS, MinCost, AvgCost, InitFunction, CostFunction, FeasibleFunction, ...
    MaxParValue, MinParValue, Population] = Init(DisplayFlag, ProblemFunction);

OPTIONS.lambda = 10; % number of offspring to produce each generation
Keep = 2; % elitism parameter: how many of the best individuals to keep from one generation to the next

% ES initialization
phiCount = 0; % number of successful mutations
sigma = 1; % standard deviation for changing solutions

% Begin the optimization loop
for GenIndex = 1 : OPTIONS.Maxgen
    % Produce children via uniform crossover
    for k = 1 : OPTIONS.lambda
        % Randomly select two parents (this works better if more fit parents are selected)
        p(1) = 1 + floor(rand * OPTIONS.popsize);
        p(2) = 1 + floor(rand * OPTIONS.popsize);
        for i = 1 : OPTIONS.numVar
            parent = 1 + (rand < 0.5);
            Population(OPTIONS.popsize+k).chrom(i) = Population(p(parent)).chrom(i);
        end
    end
    % Save the old cost values
    for i = Keep+1 : OPTIONS.popsize
        OldCost(i) = Population(i).cost;
    end
    % Mutate each individual, except don't mutate the elites.
    OPTIONS.popsize = OPTIONS.popsize + OPTIONS.lambda;
    for i = Keep+1 : OPTIONS.popsize
        Population(i).chrom = Population(i).chrom + sigma .* randn(1, OPTIONS.numVar);
    end
    % Make sure the population does not have duplicates. 
    Population = ClearDups(Population, MaxParValue, MinParValue);
    % Make sure each individual is legal.
    Population = FeasibleFunction(OPTIONS, Population);
    % Calculate cost
    Population = CostFunction(OPTIONS, Population);
    % Count how many of the mutations resulted in improvements.
    for i = Keep+1 : OPTIONS.popsize-OPTIONS.lambda
        if Population(i).cost < OldCost(i)
            phiCount = phiCount + 1;
        end
    end
    % Update sigma as required.
    MutationCount = GenIndex * (OPTIONS.popsize - OPTIONS.lambda - Keep);
    if (phiCount < MutationCount/5)
        sigma = sigma / 1.22;
    elseif (phiCount > MutationCount/5)
        sigma = sigma * 1.22;
    end
    % Sort from best to worst
    Population = PopSort(Population);
    % Eliminate the worst individuals;
    OPTIONS.popsize = OPTIONS.popsize - OPTIONS.lambda;
    Population = Population(1:OPTIONS.popsize);
    % Compute the average cost of the valid individuals
    [AverageCost, nLegal] = ComputeAveCost(Population);
    % Display info to screen
    MinCost = [MinCost Population(1).cost];
    AvgCost = [AvgCost AverageCost];
    if DisplayFlag
        disp(['The best and mean of Generation # ', num2str(GenIndex), ' are ',...
            num2str(MinCost(end)), ' and ', num2str(AvgCost(end))]);
    end
end % end of optimization loop
Conclude(DisplayFlag, OPTIONS, Population, nLegal, MinCost);
return;
