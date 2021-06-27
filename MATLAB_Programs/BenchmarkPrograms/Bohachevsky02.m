function [InitFunction, CostFunction, FeasibleFunction] = Bohachevsky02

InitFunction = @Bohachevsky02Init;
CostFunction = @Bohachevsky02Cost;
FeasibleFunction = @Bohachevsky02Feasible;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MaxParValue, MinParValue, Population, OPTIONS] = Bohachevsky02Init(OPTIONS)

global MinParValue MaxParValue
Granularity = 0.1;
MinParValue = 1;
MaxParValue = floor(1 + 2 * 10 / Granularity);
% Initialize population
for popindex = 1 : OPTIONS.popsize
    chrom = floor(MinParValue + (MaxParValue - MinParValue + 1) * rand(1,OPTIONS.numVar));
    Population(popindex).chrom = chrom;
end
OPTIONS.OrderDependent = false;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = Bohachevsky02Cost(OPTIONS, Population)

% Compute the cost of each member in Population

global MinParValue MaxParValue
popsize = OPTIONS.popsize;
p = OPTIONS.numVar;
for popindex = 1 : popsize
    sum1 = 0;
    sum2 = 0;
    Population(popindex).cost=0;
    for i = 1 : p
        gene = Population(popindex).chrom(i);
        x = (gene - MinParValue) / (MaxParValue - MinParValue) * 2 * 10 - 10;
        sum1 = sum1+(x^2)/4000;
        sum2 = sum2*cos(x/(sqrt(i)));
        Population(popindex).cost = Population(popindex).cost+sum1-sum2+1;    
    end 
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = Bohachevsky02Feasible(OPTIONS, Population)

global MinParValue MaxParValue
for i = 1 : OPTIONS.popsize
    for k = 1 : OPTIONS.numVar
        Population(i).chrom(k) = max(Population(i).chrom(k), MinParValue);
        Population(i).chrom(k) = min(Population(i).chrom(k), MaxParValue);
    end
end
return;