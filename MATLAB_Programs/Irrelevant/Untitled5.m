clear;
clc;

OptimizationFunction = char({
    %     'ACO', % ant colony optimization
    'BBO FunctionOptionsExpressions(1,:,1,3)',% biogeography-based optimization
    'BBO FunctionOptionsExpressions(1,:,2,3)',% biogeography-based optimization
    %     'DE', % differential evolution
    %     'ES', % evolutionary strategy
    %     'GA', % genetic algorithm
    %     'PBIL', % probability based incremental learning
    %     'PSO', % particle swarm optimization
    %     'StudGA', % stud genetic algorithm
    });

for i=1:size(OptimizationFunction,1)
    string010(i,:)=strsplit(OptimizationFunction(i,:));
    string020(i,:)=char(string010(i,end));
end