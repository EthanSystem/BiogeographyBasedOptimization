function [MinCost] = ACO(ProblemFunction, DisplayFlag)

% Ant colony optimization algorithm for optimizing a general function.

% INPUTS: ProblemFunction is the handle of the function that returns 
%         the handles of the initialization, cost, and feasibility functions.
%         DisplayFlag says whether or not to display information during iterations and plot results.

if ~exist('DisplayFlag', 'var')
    DisplayFlag = true;
end

[OPTIONS, MinCost, AvgCost, InitFunction, CostFunction, FeasibleFunction, ...
    MaxParValue, MinParValue, Population] = Init(DisplayFlag, ProblemFunction);

Keep = 2; % elitism parameter: how many of the best individuals to keep from one generation to the next

% ACO parameter initialization
tau0 = 1e-6; % initial pheromone value, between 0 and 0.5
Q = 20; % pheromonone update constant, between 0 and 100
q0 = 1; % exploration constant, between 0 and 1
rhog = 0.9; % global pheromone decay rate, between 0 and 1
rhol = 0.5; % local pheromone decay rate, between 0 and 1
alpha = 1; % pheromone sensitivity, between 1 and 5
beta = 5; % visibility sensitivity, between 0 and 15
tau = tau0 * ones(MaxParValue-MinParValue+1, 1); % initial pheromone values
p = zeros(size(tau)); % allocate array for probabilities

% Begin the optimization loop
for GenIndex = 1 : OPTIONS.Maxgen
    % pheromone decay
    tau = (1 - rhog) * tau;
    % Use each solution to update the pheromone for each parameter value
    for k = 1 : OPTIONS.popsize
        Cost = Population(k).cost;
        Chrom = Population(k).chrom;
        for i = 1 : length(Chrom)
            j = Chrom(i);
            if (Cost == 0)
                tau(j-MinParValue+1) = max(tau);
            else
                tau(j-MinParValue+1) = tau(j-MinParValue+1) + Q / Cost;
            end
        end    
    end
    % Use the probabilities to generate new solutions
    for k = Keep+1 : OPTIONS.popsize
        for j = 1 : OPTIONS.numVar
            % Generate probabilities based on pheromone amounts
            p = tau .^ alpha;
            p = p / sum(p);
            [Maxp, Maxpindex] = max(p);
            if rand < q0
                Select_index = Maxpindex;
            else
                SelectProb = p(1);
                Select_index = 1;
                RandomNumber = rand;
                while SelectProb < RandomNumber
                    Select_index = Select_index + 1;
                    if Select_index >= MaxParValue - MinParValue + 1
                        break;
                    end
                    SelectProb = SelectProb + p(Select_index);
                end
            end
            Population(k).chrom(j) = MinParValue + Select_index - 1;
            % local pheromone update
            tau(Select_index) = (1 - rhol) * tau(Select_index) + rhol * tau0;     
        end
    end
    % Make sure the population does not have duplicates. 
    Population = ClearDups(Population, MaxParValue, MinParValue);
    % Make sure each individual is legal.
    Population = FeasibleFunction(OPTIONS, Population);
    % Calculate cost
    Population = CostFunction(OPTIONS, Population);
    % Sort from best to worst
    Population = PopSort(Population);
    % Compute the average cost of the valid individuals
    [AverageCost, nLegal] = ComputeAveCost(Population);
    % Display info to screen
    MinCost = [MinCost Population(1).cost];
    AvgCost = [AvgCost AverageCost];
    if DisplayFlag
        disp(['The best and mean of Generation # ', num2str(GenIndex), ' are ',...
            num2str(MinCost(end)), ' and ', num2str(AvgCost(end))]);
    end
end
Conclude(DisplayFlag, OPTIONS, Population, nLegal, MinCost);
return;
