function [InitFunction, CostFunction, FeasibleFunction] = HartmanFamily19

InitFunction = @HartmanFamily19Init;
CostFunction = @HartmanFamily19Cost;
FeasibleFunction = @HartmanFamily19Feasible;
%原则上只考虑n=3个变量数
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MaxParValue, MinParValue, Population, OPTIONS] = HartmanFamily19Init(OPTIONS)

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
function [Population] = HartmanFamily19Cost(OPTIONS, Population)

% Compute the cost of each member in Population
%原则上只考虑n=3个变量数
global MinParValue MaxParValue
popsize = OPTIONS.popsize;
p = OPTIONS.numVar;
m=4;n=3;%函数中参数设置
a=[3 10 30;0.1 10 35;3 10 30;0.1 10 35];
c=[1 1.2 3 3.2]';
pp=[0.3689 0.1170 0.2673;0.4699 0.4387 0.7470;
   0.1091 0.8732 0.5547;0.038150 0.5743 0.8828];
for popindex = 1 : popsize
    Population(popindex).cost=0;
    for i = 1 : p
        gene = Population(popindex).chrom(i);
        x(i) = (gene - MinParValue) / (MaxParValue - MinParValue) * 2 * 0.5 - 0.5;
    end
    sum2=0;
    for i=1:m
        sum1=0;
        for j=1:n
            sum1=sum1-a(i,j)*((x(j)-pp(i,j))^2);
        end
%         sum2=sum2-c(i)*exp(sum1);
    end
    Population(popindex).cost = sum2;    
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = HartmanFamily19Feasible(OPTIONS, Population)

global MinParValue MaxParValue
for i = 1 : OPTIONS.popsize
    for k = 1 : OPTIONS.numVar
        Population(i).chrom(k) = max(Population(i).chrom(k), MinParValue);
        Population(i).chrom(k) = min(Population(i).chrom(k), MaxParValue);
    end
end
return;