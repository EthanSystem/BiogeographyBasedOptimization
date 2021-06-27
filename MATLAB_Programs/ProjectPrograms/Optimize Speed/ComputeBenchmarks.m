function [ Population,cost ] = ComputeBenchmarks( ProblemFunction,OPTIONS )
%COMPUTEBENCHMARKS Summary of this function goes here
%   Detailed explanation goes here



% Get the addresses of the initialization, cost, and feasibility functions.
[InitFunction, CostFunction, FeasibleFunction] = ProblemFunction();
% Initialize the population.
[MaxParValue, MinParValue, Population, OPTIONS] = InitFunction(OPTIONS);
% Make sure each individual is legal.
Population = FeasibleFunction(OPTIONS, Population);
% Calculate cost
Population = CostFunction(OPTIONS, Population);

cost=cat(1,Population.cost);

return;
end

