function [InitFunction, CostFunction, FeasibleFunction] = BraninRCOS
InitFunction = @BraninRCOSInit;
CostFunction = @BraninRCOSCost;
FeasibleFunction = @BraninRCOSFeasible;
%原则上只允许一个解中只有两个变量
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MaxParValue, MinParValue, Population, OPTIONS] = BraninRCOSInit(OPTIONS)

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
function [Population] = BraninRCOSCost(OPTIONS, Population)

% Compute the cost of each member in Population
%原则上只允许一个解中只有两个变量，OPTIONS.numVar=2
global MinParValue MaxParValue
popsize = OPTIONS.popsize;
p = OPTIONS.numVar;
for popindex = 1 : popsize
    Population(popindex).cost=0;sum1=0;sum2=0;
        for i = 1 : p
            gene = Population(popindex).chrom(i);
            x(i) = (gene - MinParValue) / (MaxParValue - MinParValue) * 2 * 10 - 10;
        end
        sum1=(x(2)-(5.1/(4*pi^2)))*x(1)^2+(5/pi*x(1)-6)^2;
        sum2=10*(1-(1/(8*pi)))*cos(x(1))+10;
       Population(popindex).cost =sum1+sum2; 
end 
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = BraninRCOSFeasible(OPTIONS, Population)

global MinParValue MaxParValue
for i = 1 : OPTIONS.popsize
    for k = 1 : OPTIONS.numVar
        Population(i).chrom(k) = max(Population(i).chrom(k), MinParValue);
        Population(i).chrom(k) = min(Population(i).chrom(k), MaxParValue);
    end
end
return;