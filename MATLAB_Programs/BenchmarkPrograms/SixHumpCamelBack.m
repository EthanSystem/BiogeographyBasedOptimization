function [InitFunction, CostFunction, FeasibleFunction] = SixHumpCamelBack

InitFunction = @SixHumpCamelBackInit;
CostFunction = @SixHumpCamelBackCost;
FeasibleFunction = @SixHumpCamelBackFeasible;
%原则上只允许一个解中只有两个变量
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MaxParValue, MinParValue, Population, OPTIONS] = SixHumpCamelBackInit(OPTIONS)

global MinParValue MaxParValue
Granularity = 0.01;
MinParValue = 1;
MaxParValue = floor(1 + 2 * 5 / Granularity);
% Initialize population
for popindex = 1 : OPTIONS.popsize
    chrom = floor(MinParValue + (MaxParValue - MinParValue + 1) * rand(1,OPTIONS.numVar));
    Population(popindex).chrom = chrom;
end
OPTIONS.OrderDependent = false;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = SixHumpCamelBackCost(OPTIONS, Population)

% Compute the cost of each member in Population
%原则上只允许一个解中只有两个变量，OPTIONS.numVar=2
global MinParValue MaxParValue
popsize = OPTIONS.popsize;
p = OPTIONS.numVar;
for popindex = 1 : popsize
    Population(popindex).cost=0;
        for i = 1 : p
            gene = Population(popindex).chrom(i);
            x(i) = (gene - MinParValue) / (MaxParValue - MinParValue) * 2 * 5 - 5;
        end
        Population(popindex).cost =4*x(1)^2-2.1*x(1)^4+(1/3)*x(1)^6+x(1)*x(2)-4*x(2)^2+4*x(2)^4;
      
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = SixHumpCamelBackFeasible(OPTIONS, Population)

global MinParValue MaxParValue
for i = 1 : OPTIONS.popsize
    for k = 1 : OPTIONS.numVar
        Population(i).chrom(k) = max(Population(i).chrom(k), MinParValue);
        Population(i).chrom(k) = min(Population(i).chrom(k), MaxParValue);
    end
end
return;