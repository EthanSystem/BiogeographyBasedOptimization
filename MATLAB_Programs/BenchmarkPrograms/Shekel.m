function [InitFunction, CostFunction, FeasibleFunction] = Shekel
InitFunction = @ShekelInit;
CostFunction = @ShekelCost;
FeasibleFunction = @ShekelFeasible;
%ԭ����ֻ����һ������ֻ����������
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MaxParValue, MinParValue, Population, OPTIONS] = ShekelInit(OPTIONS)

global MinParValue MaxParValue
Granularity = 0.1;
MinParValue = 1;
MaxParValue = floor(1 + 2 * 65.536 / Granularity);
% Initialize population
for popindex = 1 : OPTIONS.popsize
    chrom = floor(MinParValue + (MaxParValue - MinParValue + 1) * rand(1,OPTIONS.numVar));
    Population(popindex).chrom = chrom;
end
OPTIONS.OrderDependent = false;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = ShekelCost(OPTIONS, Population)

% Compute the cost of each member in Population
%ԭ����ֻ����һ������ֻ������������OPTIONS.numVar=2
global MinParValue MaxParValue
popsize = OPTIONS.popsize;
p = OPTIONS.numVar;
a=[];k=-32;kk=-32;
bnum=25;%�������в���ֵ25
for j=1:bnum
    if(~mod(j,5))
       a(1,j)=kk;
        kk=-32;%���Ƶ�һ��
        a(2,j)=k;
        k=k+16;%���Ƶڶ���
    else
        a(1,j)=kk;
        kk=kk+16;
        a(2,j)=k;
    end

end%�õ���������ֵa(i,j)
for popindex = 1 : popsize
    Population(popindex).cost=0;
    for i = 1 : p
        gene = Population(popindex).chrom(i);
        x(i) = (gene - MinParValue) / (MaxParValue - MinParValue) * 2 * 10 - 10;
    end
    sum2=0;
    for j=1:25
        i=1;
        sum1=(x(i)-a(i,j))^6+(x(i+1)-a(i+1,j))^6;
        sum2=sum2+(j+sum1)^-1;    
    end
    Population(popindex).cost = (1/500+sum2)^-1;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = ShekelFeasible(OPTIONS, Population)

global MinParValue MaxParValue
for i = 1 : OPTIONS.popsize
    for k = 1 : OPTIONS.numVar
        Population(i).chrom(k) = max(Population(i).chrom(k), MinParValue);
        Population(i).chrom(k) = min(Population(i).chrom(k), MaxParValue);
    end
end
return;