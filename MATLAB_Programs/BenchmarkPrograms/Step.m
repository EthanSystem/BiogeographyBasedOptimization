function [InitFunction, CostFunction, FeasibleFunction] = Step(OPTIONS)

InitFunction = @Step_ameliorateInit;
CostFunction = @Step_ameliorateCost;
FeasibleFunction = @Step_ameliorateFeasible;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MaxParValue, MinParValue, Population, OPTIONS] = Step_ameliorateInit(OPTIONS)

global MinParValue MaxParValue
Granularity = 0.1;
MinParValue = 1;
MaxParValue = floor(1 + 200 / Granularity);
% Initialize population
for popindex = 1 : OPTIONS.popsize
    %% 验证函数正确性用 %%%%%%%%
    %     chrom = floor(MinParValue + (MaxParValue - MinParValue + 1) * OPTIONS.randNum(popindex,:));
    %%%%%%%%%%%%%%%%%%%%%
    chrom = floor(MinParValue + (MaxParValue - MinParValue + 1) * rand(1,OPTIONS.numVar));
    Population(popindex).chrom = chrom;
    Population(popindex).cost = 0;
end
OPTIONS.OrderDependent = false;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = Step_ameliorateCost(OPTIONS, Population)
% Compute the cost of each member in Population
global MinParValue MaxParValue

gene=cat(1,Population(:).chrom);

x = (gene - MinParValue) / (MaxParValue - MinParValue) * 2 * 100 - 100;
t=floor(x+0.5).^2;
F=sum(t,2);


for i=1:OPTIONS.popsize
    Population(i).cost=F(i,:);
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = Step_ameliorateFeasible(OPTIONS, Population)

global MinParValue MaxParValue
for i = 1 : OPTIONS.popsize
    Population(i).chrom = max(Population(i).chrom, MinParValue);
    Population(i).chrom = min(Population(i).chrom, MaxParValue);
end
return;