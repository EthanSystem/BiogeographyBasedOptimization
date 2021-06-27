function [MinCost] = DE(ProblemFunction, DisplayFlag)

% Differential evolution algorithm for optimizing a general function.

% INPUTS: ProblemFunction is the handle of the function that returns
%         the handles of the initialization, cost, and feasibility functions.
%         DisplayFlag says whether or not to display information during iterations and plot results.

if ~exist('DisplayFlag', 'var')
    DisplayFlag = true;
end

[OPTIONS, MinCost, AvgCost, InitFunction, CostFunction, FeasibleFunction, ...
    MaxParValue, MinParValue, Population] = Init(DisplayFlag, ProblemFunction);

% DE parameter initialization
F = 0.5; % weighting factor
CR = 0.5; % crossover constant
PopTest = Population(1:2); % set aside space for temporary trial/target population structure

% Begin the optimization loop
for GenIndex = 1 : OPTIONS.Maxgen
    for k = 1 : OPTIONS.popsize
        % Generate the mutant "v"
        r1 = round(OPTIONS.popsize * rand + 0.5);
        while true
            r2 = round(OPTIONS.popsize * rand + 0.5);
            if (r2 ~= r1), break, end
        end
        while true
            r3 = round(OPTIONS.popsize * rand + 0.5);
            if (r3 ~= r1) & (r3 ~= r2), break, end
        end
        v.chrom = Population(r1).chrom + F * (Population(r2).chrom - Population(r3).chrom);
        % Select the index "r4" for the target vector "x"
        r4 = round(OPTIONS.popsize * rand + 0.5);
        % Generate the trial vector "u"
        for j = 1 : OPTIONS.numVar
            if rand < CR
                uchrom(j) = v.chrom(j);
            else
                uchrom(j) = Population(r4).chrom(j);
            end
        end
        % Create a small two-member population to decide whether to keep the trial vector or target vector
        PopTest(1).chrom = uchrom;
        PopTest(2) = Population(r4);
        SavePopSize = OPTIONS.popsize;
        OPTIONS.popsize = 2;
        % Make sure the trial vector is feasible
        PopTest = FeasibleFunction(OPTIONS, PopTest);
        % Compute the cost of the trial vector and the test vector
        PopTest = CostFunction(OPTIONS, PopTest);
        % Restore the original population size
        OPTIONS.popsize = SavePopSize;
        % Decide whether or not to replace the target vector with the trial vector
        if PopTest(1).cost < Population(r4).cost
            Population(r4) = PopTest(1);
        end
    end
    % Make sure the population does not have duplicates.
    Population = ClearDups(Population, MaxParValue, MinParValue);
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
