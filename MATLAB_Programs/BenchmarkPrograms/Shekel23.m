function [InitFunction, CostFunction, FeasibleFunction] = Shekel23
InitFunction = @Shekel23Init;
CostFunction = @Shekel23Cost;
FeasibleFunction = @Shekel23Feasible;
%原则上只允许一个解中只有10个变量
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MaxParValue, MinParValue, Population, OPTIONS] = Shekel23Init(OPTIONS)

global MinParValue MaxParValue
Granularity = 0.1;
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
function [Population] = Shekel23Cost(OPTIONS, Population)

% Compute the cost of each member in Population
%原则上只允许一个解中只有10个变量，
global MinParValue MaxParValue
popsize = OPTIONS.popsize;
p = 10;%OPTIONS.numVar=5
a=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;
    3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;
    6 2 6 2;7 3.6 7 3.6];
c=[0.1 0.2 0.2 0.4 0.4 0.6 0.3 0.7 0.5 0.5]';
for popindex = 1 : popsize
    Population(popindex).cost=0;
    sum=0;
    for i = 1 : p
        gene = Population(popindex).chrom(i);
        x = (gene - MinParValue) / (MaxParValue - MinParValue) * 2 * 5 - 5;
        sum=sum-((x-a(i)*(x-a(i))+c(i)));
    end
    Population(popindex).cost =sum ;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = Shekel23Feasible(OPTIONS, Population)

global MinParValue MaxParValue
for i = 1 : OPTIONS.popsize
    for k = 1 : OPTIONS.numVar
        Population(i).chrom(k) = max(Population(i).chrom(k), MinParValue);
        Population(i).chrom(k) = min(Population(i).chrom(k), MaxParValue);
    end
end
return;