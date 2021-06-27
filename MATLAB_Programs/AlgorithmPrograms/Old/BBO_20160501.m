function [MinCost, Hamming, ExportData] = BBO(SET,OPTIONS, ExportData, ProblemFunction , ComputeFlagOfFitnessFunction ,  OperatorFlagOfMigration ,ComputeFlagOfMigration, ComputeFlagOfMutation , DisplayFlag , RandSeed)
% Biogeography-based optimization (BBO) software for minimizing a general function
% BBO算法用于最小化通用函数的解的误差值
% INPUTS: ProblemFunction is the handle of the function that returns
%         the handles of the initialization, cost, and feasibility
%         functions.  初始化、价值、可行性函数的句柄
%        ComputeFlagOfFitnessFunction = { 'species count fitness',
%        'normalized fitness' } . 分别表示： '使用层次适应值（也就是物种数）作为适应值，这个是原作者的思路' , '使用通过函数值直接获得的，归一化的适应值' 。
%         OperatorFlagOfMigration = { 'PI', 'TI', 'PE', 'TE' } . They means partial immigration-based BBO , total immigration-based BBO , partial emigration-based BBO , total emigration-based BBO respectively .
%         ComputeFlagOfMigration = { 'species count fitness' , 'normalized fitness' } . They
%           means ：' 使用层次适应值（也就是物种数）作为适应值，这个是原作者的思路 ' , '
%           使用通过函数值直接获得的，归一化的适应值 ' 。
%          ComputeFlagOfMutation = { 'by prob' , 'by steady' , 'by constant' } . 分别表示： '当物种数分布概率与物种数的变化有关时，计算种群的变异率' ,
%          '当物种数分布概率仅与迁入迁出率有关而与物种数的变化无关时，计算种群的变异率' ,
%          '当物种数分布概率呈固定值时，计算种群的变异率' 。
%         DisplayFlag = true or false, whether or not to display and plot results.
%         ProbFlag = true or false, whether or not to use probabilities to update emigration rates.
%         RandSeed = random number seed
% OUTPUTS: MinCost = array of best solution, one element for each generation
%          Hamming = final Hamming distance between solutions
%          解与解之间的Hamming距离
% CAVEAT: The "ClearDups" function that is called below replaces duplicates with randomly-generated
%         individuals, but it does not then recalculate the cost of the replaced individuals.


if ~exist('DisplayFlag', 'var')
    DisplayFlag = true;
end
if ~exist('ComputeFlagOfFitnessFunction', 'var')
    ComputeFlagOfFitnessFunction = 'species count fitness';
end
if ~exist('OperatorFlagOfMigration', 'var')
    OperatorFlagOfMigration = 'PI';
end

if ~exist('ComputeFlagOfMigration', 'var')
    ComputeFlagOfMigration = 'L-L';
end
if ~exist('ComputeFlagOfMutation', 'var')
    ComputeFlagOfMutation = 'by prob';
end
if ~exist('RandSeed', 'var')
    RandSeed = round(sum(100*clock));
end

[OPTIONS, MinCost, AvgCost, InitFunction, CostFunction, FeasibleFunction, ...
    MaxParValue, MinParValue, Population, ExportData] = Init(SET, OPTIONS, DisplayFlag, ExportData, ProblemFunction, RandSeed);

Population = CostFunction(OPTIONS, Population);

OPTIONS.pmodify = 1; % habitat modification probability 栖息地的迁移概率

OPTIONS.pmutate = 0.001; % initial mutation probability初始化突变概率

Keep = 2; % elitism parameter: how many of the best habitats to keep from one generation to the next 精英参数
lambdaLower = 0.0; % lower bound for immigration probabilty per gene
lambdaUpper = 1.0; % upper bound for immigration probabilty per gene
muLower = 0.0; % lower bound for emigration probability per gene
muUpper = 1.0; % upper bound for emigration probablity per gene
dt = 1; % step size used for numerical integration of probabilities
I = 1; % max immigration rate for each island
E = 1; % max emigration rate, for each island
P = OPTIONS.popsize; % max species count, for each island



% Initialize the species count probability of each habitat
% Later we might want to initialize probabilities based on cost
%   初始化每一个栖息地的种群数量的概率
for j = 1 : length(Population)
    Prob(j) = 1 / length(Population);
end




%% Begin the optimization loop
for GenIndex = 1 : OPTIONS.Maxgen   % 从第一代开始
    %% Save the best habitats in a temporary array.
    for j = 1 : Keep
        chromKeep(j,:) = Population(j).chrom;
        costKeep(j) = Population(j).cost;
    end
    
    %% 计算适应值函数 Map cost values to species counts.
    [Population] = GetFitnessCost(ComputeFlagOfFitnessFunction,Population, P);
    
    %% Migration 选取迁移率模型，计算迁移率
    % Compute immigration rate and emigration rate for each species count.
    % lambda(i) is the immigration rate for habitat i.
    % mu(i) is the emigration rate for habitat i.
    [lambda, mu] = GetLambdaMu(ComputeFlagOfMigration,Population, I, E, P);
    
    
    %% 设置两种不同的种群分布概率计算方式
    % Compute the time derivative of Prob(i) for each habitat i.
    % 计算每个栖息地的概率对时间的导数
    % 注意：再次强调，运行函数前，已经假定栖息地 i已经从优到劣排序好了。
    switch ComputeFlagOfMutation         % 依据层次适应值k计算迁出率 If allowed to use probablities to update emigration rates .
        case 'by prob'         % 计算每个栖息地的概率对时间的导数
            Prob=GetProbablisticProb(Population,Prob,I,E,P,lambda,mu,dt);
        case 'by steady'     % 迁移率设为常数
            Prob=GetSteadyStateProb (Population,Prob,I,lambda,mu);
    end
    
    %% 设置 BBO 迁移方式的算子：
    [Island]=MigrationMethod(OPTIONS, OperatorFlagOfMigration, Population, lambda, mu, lambdaLower, lambdaUpper, muLower, muUpper);
    
    %%  Mutation 计算变异率
    [Island]=ComputeMutation(OPTIONS, ComputeFlagOfMutation, Population, Island, Prob,MinParValue, MaxParValue);
    
    %% 更新
    % Replace the habitats with their new versions.
    for k = 1 : length(Population)
        Population(k).chrom = Island(k,:);
    end
    % Make sure each individual is legal.
    Population = FeasibleFunction(OPTIONS, Population);
    % Calculate cost
    Population = CostFunction(OPTIONS, Population);
    % Sort from best to worst
    Population = PopSort(Population);
    % Replace the worst with the previous generation's elites.
    n = length(Population);
    for k = 1 : Keep
        Population(n-k+1).chrom = chromKeep(k,:);
        Population(n-k+1).cost = costKeep(k);
    end
    % Make sure the population does not have duplicates.
    Population = ClearDups(Population, MaxParValue, MinParValue);
    % Sort from best to worst
    Population = PopSort(Population);
    % Compute the average cost
    [AverageCost, nLegal] = ComputeAveCost(Population);
    % Display info to screen
    MinCost = [MinCost Population(1).cost];
    AvgCost = [AvgCost AverageCost];
    
    %% 显示与输出
    if DisplayFlag
        disp(['The best and mean of Generation # ', num2str(GenIndex), ' are ',...
            num2str(MinCost(end)), ' and ', num2str(AvgCost(end))]);
    end
    ExportData.BestAndMeanOfGeneration(GenIndex+1).generation(SET.monte_index)=GenIndex;
    ExportData.BestAndMeanOfGeneration(GenIndex+1).best(SET.monte_index)=MinCost(end);
    ExportData.BestAndMeanOfGeneration(GenIndex+1).mean(SET.monte_index)=AvgCost(end);

end

%% 结论
[ExportData]=Conclude(SET, OPTIONS, DisplayFlag, ExportData, Population, nLegal, MinCost);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%可以改进
%% Plot some results
% PlotConcludeFigures( SET, OPTIONS , 'BBO' , MinCost)

    
%% Obtain a measure of population diversity
% 从栖息地k到下一层级基因的遍历
for k = 1 : length(Population)
    Chrom = Population(k).chrom;
    for j = MinParValue : MaxParValue
        indices = find(Chrom == j);
        CountArr(k,j) = length(indices); % array containing gene counts of each habitat
    end
end
Hamming = 0;
for m = 1 : length(Population)
    for j = m+1 : length(Population)
        for k = MinParValue : MaxParValue
            Hamming = Hamming + abs(CountArr(m,k) - CountArr(j,k));
        end
    end
end

if DisplayFlag
    disp(['Diversity measure = ', num2str(Hamming)]);
end
ExportData.DiversityMeasure(SET.monte_index)=Hamming;

return;



end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GetFitnessCost 函数
function [Population] = GetFitnessCost(ComputeFlagOfFitnessFunction,Population, P)
% Map cost values to fitness values.
% This loop assumes the population is already sorted from most fit to least fit.
% 循环前，个体的值已经从好到坏排好序
switch ComputeFlagOfFitnessFunction
    % Map cost values to species counts.
    % 使用层次适应值（也就是物种数）作为适应值，这个是原作者的思路
    case 'species count fitness'
        for i = 1 : length(Population)
            if Population(i).cost < inf
                Population(i).FitnessCost = P - i;
            else
                Population(i).FitnessCost = 0;
            end
        end
        
        % 使用通过函数值直接获得的，归一化的适应值
    case 'normalized fitness'
        minCost=Population(1).cost;
        maxCost=Population(length(Population)).cost;
        %         temp=ones(length(Population),1);
        for i = 1 : length(Population)
            temp(i)=-Population(i).cost;
        end
        for i = 1 : length(Population)
            if Population(i).cost < inf
                Population(i).FitnessCost = (temp(i)+maxCost)./(maxCost-minCost);
            else
                Population(i).FitnessCost = 0;
            end
        end
        
end
return;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GetLambdaMu函数
function [lambda, mu] = GetLambdaMu(ComputeFlagOfMigration,Population, I, E, P)
switch ComputeFlagOfMigration
    case 'L-L'
        % 迁移模型采用 线性迁入-线性迁出 方式
        % Compute immigration rate and extinction rate for each species count.
        % lambda(i) is the immigration rate for individual i.
        % mu(i) is the extinction rate for individual i.
        
        for i = 1 : length(Population)
            lambda(i) = I * (1 - Population(i).FitnessCost / P);
            mu(i) = E * Population(i).FitnessCost / P;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%% 待以后扩充
        
end

return;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GetSteadyStateOfProb 函数
function [Prob]=GetSteadyStateProb (Population,Prob,I,lambda,mu)
% 仅与迁入迁出率有关而与物种数的变化无关时的物种数分布概率
% calculate steadyProb(0)
q=0;
p=zeros(1,length(Population));
p(1)=I/mu(1);
for i=2:length(Population)
    p(i)=p(i-1)*lambda(i-1)/mu(i);
end
for i=1:length(Population)
    q=q+p(i);
end
Prob0=1/(1+q);
% calculate steadyProb(i)
for i=1:length(Population)
    Prob(i)=Prob0*p(i);
end

return;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GetProbablisticProb 函数
function [Prob]=GetProbablisticProb(Population,Prob,I,E,P,lambda,mu,dt)
%% Compute the time derivative of Prob(i) for each habitat i.
% 计算每个栖息地的概率对时间的导数
% 注意：再次强调，运行函数前，已经假定栖息地 i已经从优到劣排序好了。
% ProbDot=zeros(1:length(Population));
for j = 1 : length(Population)
    % Compute lambda for one less than the species count of habitat i.
    lambdaMinus = I * (1 - (Population(j).FitnessCost - 1) / P);
    % Compute mu for one more than the species count of habitat i.
    muPlus = E * (Population(j).FitnessCost + 1) / P;
    % Compute Prob for one less than and one more than the species count of habitat i.
    % Note that species counts are arranged in an order opposite to that presented in
    % MacArthur and Wilson's book - that is, the most fit
    % habitat has index 1, which has the highest species count.
    if j < length(Population)
        ProbMinus = Prob(j+1);
    else
        ProbMinus = 0;
    end
    if j > 1
        ProbPlus = Prob(j-1);
    else
        ProbPlus = 0;
    end
    ProbDot(j) = -(lambda(j) + mu(j)) * Prob(j) + lambdaMinus * ProbMinus + muPlus * ProbPlus;
end

% Compute the new probabilities for each species count.
Prob = Prob + ProbDot * dt;
Prob = max(Prob, 0);    % 使Prob的所有值非负
Prob = Prob / sum(Prob);
return;
end


%% MigrationMethod 迁移算子函数
function [Island]=MigrationMethod(OPTIONS, OperatorFlagOfMigration, Population, lambda, mu, lambdaLower, lambdaUpper, muLower, muUpper)
%% 设置 BBO 迁移方式的算子：
switch OperatorFlagOfMigration
    case 'PI'  % it means partial immigration-based BBO       
        %% Now use lambda and mu to decide how much information to share between habitats.
        % 迁移运算：迁入主导部分变量迁移。用 lambda 和 mu 确定在 habitats 之间可以分享的信息量。
        lambdaMin = min(lambda);
        lambdaMax = max(lambda);
        for k = 1 : length(Population)
            if rand > OPTIONS.pmodify
                continue;
            end
            % Normalize the immigration rate.
            % lamdaScale 将现有的 lambda(k) 的尺度压缩到 [0,1] 实数空间，使之可以与rand()比较
            lambdaScale = lambdaLower + (lambdaUpper - lambdaLower) * (lambda(k) - lambdaMin) ./ (lambdaMax - lambdaMin);
            %% Probabilistically input new information into habitat i
            for j = 1 : OPTIONS.numVar
                if rand < lambdaScale    % 根据lambdaScale确定这个特征值是否迁徙
                    % Pick a habitat from which to obtain a feature
                    % 根据mu 选取一个 habitat 作为迁出者，通过轮盘赌选取法，mu越大的habitat越容易被选中
                    RandomNum = rand * sum(mu);
                    Select = mu(1);
                    SelectIndex = 1;
                    while (Select < RandomNum) & (SelectIndex < OPTIONS.popsize)
                        SelectIndex = SelectIndex + 1;
                        Select = Select + mu(SelectIndex);
                    end
                    %发生迁徙，接下来将把选中的栖息地的第 j 个SIV复制给第k个栖息地
                    Island(k,j) = Population(SelectIndex).chrom(j);
                else
                    Island(k,j) = Population(k).chrom(j);       % 不发生迁徙
                end
            end
        end
        

        
        
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %   算法二：矢量化算法，实测程序运行速度变慢。
%         %         输入：
%         %             Population
%         %         过程：
%         %             用于判断OPTIONS.pmodify的rand列向量randPmodify
%         %             把OPTIONS.pmodify做成列向量形式，以利于判断lambda是否迁徙的rand列向量
%         %             Population赋给In
%         %             用isPmodify对每个个体判断是否迁移，生成迁移逻辑矩阵L01，用L01更新In产生Temp01。
%         %             生成标准化迁移率矩阵lambdaScale，行=迁入个体，列=特征。用lambdaScale取代lambda
%         %             产生用于判断lambda是否迁徙的rand行向量randMigration
%         %             用lambdaScale判断该特征值是否迁移，产生迁移逻辑矩阵L02，并用L02更新In产生Temp02。
%         %             轮盘赌选择一个迁出个体，记入Index：
%         %                 建立轮盘赌指针RandomPointer
%         %                 下三角元矩阵→mu的累计矩阵Mu→求和→依次减去RandomPointer→按照正数或零变成1，负数变成0，转换成逻辑数组L03→求和+1，得到迁出个体下标Index。
%         %             迁出个体记入迁移矩阵Index：  行=迁入个体，列=特征。
%         %             用L01、L02更新Index
%         %           用上述轮盘赌过程产生Index矩阵
%         %            将Index映射到迁移关联矩阵Migration：行=迁入个体，列=迁出个体，页=特征。
%         %           用Migration实现迁出操作，产生迁出矩阵Temp03。
%         %           将迁出矩阵Temp03赋值生成 Island
%         %         输出：
%         %             Island
% 
%         lambdaMin = min(lambda);
%         lambdaMax = max(lambda);
%         randPmodify=rand(OPTIONS.popsize,1);
%         In.chrom=cat(1,Population.chrom);
%         In.cost=cat(1,Population.cost);
%         L01(randPmodify<OPTIONS.pmodify)=1;
%         L01=L01';
%         L01=repmat(L01,1,OPTIONS.numVar);
%         Temp01=In.chrom.*L01;
%         
%         lambdaScale=lambdaLower + (lambdaUpper - lambdaLower) .* (lambda - lambdaMin) ./ (lambdaMax - lambdaMin);
%         lambdaScale=lambdaScale';
%         lambdaScale=repmat(lambdaScale,1,OPTIONS.numVar);
%         randMigration=rand(OPTIONS.popsize,OPTIONS.numVar);
%         L02(randMigration<lambdaScale)=1;
%         L02=reshape(L02,OPTIONS.popsize,OPTIONS.numVar);
%         Temp02=Temp01.*L02;
%         % 轮盘赌选择一个迁出个体
%         Mu=tril(repmat(mu,length(mu),1));
%         Mu=sum(Mu,2);
%         for i=1:OPTIONS.popsize
%             for j=1:OPTIONS.numVar
%                 RandomPointer=repmat(rand*sum(mu), OPTIONS.popsize,1);
%                 distance=RandomPointer-Mu;
%                 distance(distance>=0)=1;
%                 distance(distance<0)=0;
%                 Index(i,j)=sum(distance)+1;
%             end
%         end
%         Index=Index.*L02.*L01;
%         % 产生Index矩阵
%         Migration=zeros(OPTIONS.popsize,OPTIONS.popsize,OPTIONS.numVar);
%         for i=1:OPTIONS.numVar
%             for j=1:OPTIONS.popsize
%                 if Index(j,i)==0
%                     continue;
%                 end
%                 Migration(j,Index(j,i),i)=1;
%             end
%         end
%         % 赋值过程
%         for i=1:OPTIONS.numVar
%             Temp03(:,i)=Migration(:,:,i)*Temp02(:,i);
%         end
%             Island=In.chrom;
%             Island=Island.*~Temp03;
%             Island=Island+Temp03;
        

            
            
    case 'TI' % it means total immigration-based BBO
        %% Now use lambda and mu to decide how much information to share between habitats.
        % 迁移运算：迁入主导全部变量迁移。用 lambda 和 mu 确定在 habitats 之间可以分享的信息量。
        lambdaMin = min(lambda);
        lambdaMax = max(lambda);
        for k = 1 : length(Population)
            if rand > OPTIONS.pmodify
                continue;
            end
            % Normalize the immigration rate.
            % lamdaScale 将现有的 lambda(k) 的尺度压缩到 [0,1] 实数空间，使之可以与rand()比较
            lambdaScale = lambdaLower + (lambdaUpper - lambdaLower) * (lambda(k) - lambdaMin) / (lambdaMax - lambdaMin);
            %% Probabilistically input new information into habitat i
            if rand < lambdaScale    % 根据lambdaScale确定这个特征值是否迁徙
                for j=1 : OPTIONS.numVar   % 对于每个特征值，都要选择各自的迁出者对应的特征值进行迁移
                    % 根据mu 选取一个 habitat 作为迁出者，通过轮盘赌选取法，mu越大的habitat越容易被选中
                    RandomNum = rand * sum(mu);
                    Select = mu(1);
                    SelectIndex = 1;
                    while (Select < RandomNum) & (SelectIndex < OPTIONS.popsize)
                        SelectIndex = SelectIndex + 1;
                        Select = Select + mu(SelectIndex);
                    end
                    %% 对于所有的 特征值 都发生迁徙
                    %发生迁徙，接下来将把选中的栖息地的第 j 个SIV复制给第k个栖息地
                    Island(k,j) = Population(SelectIndex).chrom(j);
                end
            else
                Island(k,:) = Population(k).chrom;       % 不发生迁徙
            end
        end
        
    case 'PE' % it means partial emigration-based BBO
        %% Now use lambda and mu to decide how much information to share between habitats.
        % 迁移运算：迁出主导部分变量迁移。用 lambda 和 mu 确定在 habitats 之间可以分享的信息量。
        muMin = min(mu);
        muMax = max(mu);
        for ii=1:length(Population)
            Island(ii,:) = Population(ii).chrom;
        end
        for k = 1 : length(Population)
            if rand > OPTIONS.pmodify
                continue;
            end
            % Normalize the emigration rate.
            % muScale 将现有的 lambda(k) 的尺度缩放到 [0,1] 实数空间，使之可以与rand()比较
            muScale = muLower + (muUpper - muLower) * (mu(k) - muMin) / (muMax - muMin);
            %% Probabilistically input new information into habitat i
            for j = 1 : OPTIONS.numVar
                if rand < muScale    % 根据lambdaScale确定这个特征值是否迁徙
                    % 根据 lambda 选取一个 habitat 作为迁入者，通过轮盘赌选取法，lambda 越大的 habitat越容易被选中
                    RandomNum = rand * sum(lambda);
                    Select = lambda(1);
                    SelectIndex = 1;
                    while (Select < RandomNum) & (SelectIndex < OPTIONS.popsize)
                        SelectIndex = SelectIndex + 1;
                        Select = Select + lambda(SelectIndex);
                    end
                    %发生迁徙，接下来将把第k个栖息地的第 j 个SIV复制给第选中的个栖息地
                    Island(SelectIndex,j) = Population(k).chrom(j);
                else
                    Island(k,j) = Population(k).chrom(j);       % 不发生迁徙
                end
            end
        end
        
    case 'TE' % it means total emigration-based BBO
        %% Now use lambda and mu to decide how much information to share between habitats.
        % 迁移运算：迁出主导全部变量迁移。用 lambda 和 mu 确定在 habitats 之间可以分享的信息量。
        muMin=min(mu);
        muMax = max(mu);
        for ii=1:length(Population)
            Island(ii,:) = Population(ii).chrom;
        end
        for k=1:length(Population)
            if rand > OPTIONS.pmodify
                continue;
            end
            % 用mu选择一个个体作为迁出者：
            muScale = muLower+(muUpper-muLower)*(mu(k)-muMin)/(muMax-muMin);        % muScale 缩放所有个体的现有的 mu 尺度到 [lower,upper] 空间，使之可以与rand()比较
            if rand<muScale
                % 对接下来选中的栖息地的所有特征值进行迁移操作
                for j=1 : OPTIONS.numVar
                    % 对该特征值，用lambda选择一个个体作为迁入者。用轮盘赌的方法。lambda 越大越容易被选中。
                    Select = lambda(1);
                    SelectIndex = 1;
                    RandomNum = rand * sum(lambda);
                    while (SelectIndex < OPTIONS.popsize ) & ( Select < RandomNum)
                        SelectIndex = SelectIndex+1;
                        Select = Select + lambda(SelectIndex);
                    end
                    % 发生迁徙，接下来把第k个栖息地的第j个SIV复制给选中的栖息地
                    Island(SelectIndex,j) = Population(k).chrom(j);
                end
            else
                Island(k,:) = Population(k).chrom;      % 不发生迁徙
            end
        end
        
end
return;
end

%%  ComputeMutation 计算变异率函数
function [Island]=ComputeMutation(OPTIONS, ComputeFlagOfMutation, Population, Island, Prob, MinParValue, MaxParValue)
switch ComputeFlagOfMutation
    case 'by prob'
        % 当物种数分布概率与物种数的变化有关时，计算种群的变异率
        Pmax = max(Prob);
        MutationRate = OPTIONS.pmutate * (1 - Prob / Pmax);
        
    case 'by steady'
        % 当物种数分布概率仅与迁入迁出率有关而与物种数的变化无关时，计算种群的变异率
        Pmax = max(Prob);
        MutationRate = OPTIONS.pmutate * (1 - Prob / Pmax);
        
    case 'by constant'
        % 当物种数分布概率呈固定值时，计算种群的变异率
        num=ones(1,OPTIONS.popsize)*OPTIONS.popsize;
        MutationRate=OPTIONS.pmutate * num;
        
end
% Mutate only the worst half of the solutions 仅让表现最差的那一半个体群发生突变
Population = PopSort(Population);
for index = round(length(Population)/2) : length(Population)    % 此时Population已经排好序了
    for parnum = 1 : OPTIONS.numVar
        if rand < MutationRate(index)
            % 新的突变值的生成采用高斯分布的生成方式
            Island(index,parnum) = floor(MinParValue + (MaxParValue - MinParValue + 1) * randn);
        end
    end
end

return;
end