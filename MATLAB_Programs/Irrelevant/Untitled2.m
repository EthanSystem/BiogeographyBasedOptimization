clear;
clc;


OPTIONS.popsize = 1; % total population size 个体总数量
OPTIONS.Maxgen = 50; % generation count limit max number of gene  生成后代代数限制 
OPTIONS.numVar = 2; % number of genes in each population member   每个个体的基因数量
OPTIONS.pmutate = 0.001; % mutation probability     突变概率

[InitFunction, CostFunction, FeasibleFunction] = Penalty1
[MaxParValue, MinParValue, Population, OPTIONS] = InitFunction(OPTIONS)
[Population] = CostFunction(OPTIONS, Population)




% Compute the cost of each member in Population

global MinParValue MaxParValue
for popindex = 1 : OPTIONS.popsize
    Population(popindex).cost = 0;
    for i = 1 : OPTIONS.numVar
        gene = Population(popindex).chrom(i);
        x = (gene - MinParValue) / (MaxParValue - MinParValue) * 2 * 50 - 50;
        y(i) = 1 + (x + 1) / 4;
        %函数中的参数值
        a=10;m=4;k=100;
        if (x > a)
            u = k * (x - 10)^m;
        elseif (x < -a)
            u = k * (-x - 10)^m;
        else
            u = 0;
        end
        Population(popindex).cost = Population(popindex).cost + u;
    end
    Population(popindex).cost = Population(popindex).cost + (10 * (sin(pi*y(1)))^2 + (y(OPTIONS.numVar) - 1)^2) * pi / OPTIONS.numVar;
    for i = 1 : OPTIONS.numVar-1
        Population(popindex).cost = Population(popindex).cost + ((y(i) - 1)^2 * (1 + 10 * (sin(pi*y(i+1)))^2)) * pi / OPTIONS.numVar;
    end
end

x=-50:1:50;
y=x;
[X,Y]=surf(x,y)









