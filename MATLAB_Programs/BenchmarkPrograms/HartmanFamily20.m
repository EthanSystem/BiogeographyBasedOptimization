function [InitFunction, CostFunction, FeasibleFunction] = HartmanFamily20

InitFunction = @HartmanFamily20Init;
CostFunction = @HartmanFamily20Cost;
FeasibleFunction = @HartmanFamily20Feasible;
%原则上每个解只考虑n=6个变量数
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MaxParValue, MinParValue, Population, OPTIONS] = HartmanFamily20Init(OPTIONS)

global MinParValue MaxParValue
Granularity = 0.01;
MinParValue = 1;
MaxParValue = floor(1 + 2 * 0.5 / Granularity);
% Initialize population
for popindex = 1 : OPTIONS.popsize
    chrom = floor(MinParValue + (MaxParValue - MinParValue + 1) * rand(1,OPTIONS.numVar));
    Population(popindex).chrom = chrom;
end
OPTIONS.OrderDependent = false;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = HartmanFamily20Cost(OPTIONS, Population)

% Compute the cost of each member in Population
%原则上每个解只考虑n=6个变量数
global MinParValue MaxParValue
popsize = OPTIONS.popsize;
p = OPTIONS.numVar;
m=4;n=6;%函数中参数设置
a=[10 3 17 3.5 1.7 8;0.05 10 17 0.1 8 14;3 3.5 1.7 10 17 8;17 8 0.05 10 0.1 14];
c=[1 1.2 3 3.2]';
pp=[0.1312 0.1696 0.5569 0.0124 0.8283 0.5886;
    0.2329 0.4135 0.8307 0.3736 0.1004 0.9991;
    0.2348 0.1415 0.3522 0.2283 0.3047 0.6650;
    0.4047 0.8828 0.8732 0.5743 0.1091 0.0381];
for popindex = 1 : popsize
    for i = 1 : p
        gene = Population(popindex).chrom(i);
        x(i) = (gene - MinParValue) / (MaxParValue - MinParValue) * 2 * 0.5 - 0.5;
    end
    sum2=0;
    for i=1:m
        sum1=0;
        for j=1:n
            sum1=sum1-a(i,j)*(x(j)-pp(i,j))^2;
        end
        sum2=sum2-c(i)*exp(sum1);
    end
    Population(popindex).cost = sum2;    
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = HartmanFamily20Feasible(OPTIONS, Population)

global MinParValue MaxParValue
for i = 1 : OPTIONS.popsize
    for k = 1 : OPTIONS.numVar
        Population(i).chrom(k) = max(Population(i).chrom(k), MinParValue);
        Population(i).chrom(k) = min(Population(i).chrom(k), MaxParValue);
    end
end
return;