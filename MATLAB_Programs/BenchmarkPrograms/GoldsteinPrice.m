function [InitFunction, CostFunction, FeasibleFunction] = GoldsteinPrice
InitFunction = @GoldsteinPriceInit;
CostFunction = @GoldsteinPriceCost;
FeasibleFunction = @GoldsteinPriceFeasible;
%原则上只允许一个解中只有两个变量
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MaxParValue, MinParValue, Population, OPTIONS] = GoldsteinPriceInit(OPTIONS)

global MinParValue MaxParValue
Granularity = 0.01;
MinParValue = 1;
MaxParValue = floor(1 + 2 * 2 / Granularity);
% Initialize population
for popindex = 1 : OPTIONS.popsize
    chrom = floor(MinParValue + (MaxParValue - MinParValue + 1) * rand(1,OPTIONS.numVar));
    Population(popindex).chrom = chrom;
end
OPTIONS.OrderDependent = false;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = GoldsteinPriceCost(OPTIONS, Population)

% Compute the cost of each member in Population
%原则上只允许一个解中只有两个变量，OPTIONS.numVar=2
global MinParValue MaxParValue
popsize = OPTIONS.popsize;
p = OPTIONS.numVar;

for popindex = 1 : popsize
    Population(popindex).cost=0;
    sumf1=0;sumf2=0;
        for i = 1 : p
            gene = Population(popindex).chrom(i);
            x(i) = (gene - MinParValue) / (MaxParValue - MinParValue) * 2 * 2 - 2;
        end
        sumf1=1+(x(1)+x(2)+1)^2*(19-14*x(1)+3*x(1)^2-14*x(2)+6*x(1)*x(2)+3*x(2)^2);
        sumf2=30+(2*x(1)-3*x(2))^2*(18-32*x(1)+12*x(1)^2+48*x(2)-36*x(1)*x(2)+27*x(2)^2);
    Population(popindex).cost =sumf1*sumf2 ;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = GoldsteinPriceFeasible(OPTIONS, Population)

global MinParValue MaxParValue
for i = 1 : OPTIONS.popsize
    for k = 1 : OPTIONS.numVar
        Population(i).chrom(k) = max(Population(i).chrom(k), MinParValue);
        Population(i).chrom(k) = min(Population(i).chrom(k), MaxParValue);
    end
end
return;