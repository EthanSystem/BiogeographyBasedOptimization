function [InitFunction, CostFunction, FeasibleFunction] = Kowalik

InitFunction = @KowalikInit;
CostFunction = @KowalikCost;
FeasibleFunction = @KowalikFeasible;
%原则上只允许一个解中只有四个变量
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MaxParValue, MinParValue, Population, OPTIONS] = KowalikInit(OPTIONS)

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
function [Population] = KowalikCost(OPTIONS, Population)

% Compute the cost of each member in Population
%原则上只允许一个解中只有四个变量，OPTIONS.numVar=2
global MinParValue MaxParValue
popsize = OPTIONS.popsize;
p = OPTIONS.numVar;
a=[0.1957 0.1947 0.1735 0.1600 0.0844 0.0627 0.0456 0.0342 0.0323 0.0235 0.0246]';
bb=[0.25 0.5 1 2 4 6 8 10 12 14 16]';
b=bb.^-1;
for popindex = 1 : popsize
    Population(popindex).cost=0;
        for i = 1 : p
            gene = Population(popindex).chrom(i);
            x(i) = (gene - MinParValue) / (MaxParValue - MinParValue) * 2 * 5 - 5;
        end
        sum1=0;sum2=0;
       for j=1:11
           sum1=x(1)*(b(j)^2+b(j)*x(2));
           sum2=b(j)^2+b(j)*x(3)+x(4);
           if(~sum2)%分为0时，不考虑后一项
                 Population(popindex).cost =Population(popindex).cost+a(j)^2;  
           else
               Population(popindex).cost =Population(popindex).cost+(a(j)-sum1/sum2)^2;
           end
       end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = KowalikFeasible(OPTIONS, Population)

global MinParValue MaxParValue
for i = 1 : OPTIONS.popsize
    for k = 1 : OPTIONS.numVar
        Population(i).chrom(k) = max(Population(i).chrom(k), MinParValue);
        Population(i).chrom(k) = min(Population(i).chrom(k), MaxParValue);
    end
end
return;