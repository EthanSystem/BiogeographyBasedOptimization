function [InitFunction, CostFunction, FeasibleFunction] = Schaffer

InitFunction = @SchafferInit;
CostFunction = @SchafferCost;
FeasibleFunction = @SchafferFeasible;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MaxParValue, MinParValue, Population, OPTIONS] = SchafferInit(OPTIONS)

global MinParValue MaxParValue
Granularity = 0.1;
MinParValue = 1;
MaxParValue = floor(1 + 2 * 3.14 / Granularity);
% Initialize population
for popindex = 1 : OPTIONS.popsize
    chrom = floor(MinParValue + (MaxParValue - MinParValue + 1) * rand(1,OPTIONS.numVar));
    Population(popindex).chrom = chrom;
end
OPTIONS.OrderDependent = false;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = SchafferCost(OPTIONS, Population)

% Compute the cost of each member in Population

global MinParValue MaxParValue
popsize = OPTIONS.popsize;
p = OPTIONS.numVar;
for popindex = 1 : popsize
    sum1 = 0;
    sum2 = 0;
    Population(popindex).cost=0;
    for i = 1 : p-1
        gene = Population(popindex).chrom(i);
        x = (gene - MinParValue) / (MaxParValue - MinParValue) * 2 * 3.14 - 3.14;
        gene2 = Population(popindex).chrom(i+1);
        y = (gene2 - MinParValue) / (MaxParValue - MinParValue) * 2 * 3.14 - 3.14;
        sum1 = sum1+x^2 + 2*y^2;
        sum2 = sum2-0.3*cos(3*pi*x)-0.4*cos(4*pi*y)+0.7;
    end
    Population(popindex).cost = sum1+sum2;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = SchafferFeasible(OPTIONS, Population)

global MinParValue MaxParValue
for i = 1 : OPTIONS.popsize
    for k = 1 : OPTIONS.numVar
        Population(i).chrom(k) = max(Population(i).chrom(k), MinParValue);
        Population(i).chrom(k) = min(Population(i).chrom(k), MaxParValue);
    end
end
return;