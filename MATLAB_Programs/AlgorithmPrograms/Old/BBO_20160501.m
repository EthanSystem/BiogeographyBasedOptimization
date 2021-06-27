function [MinCost, Hamming, ExportData] = BBO(SET,OPTIONS, ExportData, ProblemFunction , ComputeFlagOfFitnessFunction ,  OperatorFlagOfMigration ,ComputeFlagOfMigration, ComputeFlagOfMutation , DisplayFlag , RandSeed)
% Biogeography-based optimization (BBO) software for minimizing a general function
% BBO�㷨������С��ͨ�ú����Ľ�����ֵ
% INPUTS: ProblemFunction is the handle of the function that returns
%         the handles of the initialization, cost, and feasibility
%         functions.  ��ʼ������ֵ�������Ժ����ľ��
%        ComputeFlagOfFitnessFunction = { 'species count fitness',
%        'normalized fitness' } . �ֱ��ʾ�� 'ʹ�ò����Ӧֵ��Ҳ��������������Ϊ��Ӧֵ�������ԭ���ߵ�˼·' , 'ʹ��ͨ������ֱֵ�ӻ�õģ���һ������Ӧֵ' ��
%         OperatorFlagOfMigration = { 'PI', 'TI', 'PE', 'TE' } . They means partial immigration-based BBO , total immigration-based BBO , partial emigration-based BBO , total emigration-based BBO respectively .
%         ComputeFlagOfMigration = { 'species count fitness' , 'normalized fitness' } . They
%           means ��' ʹ�ò����Ӧֵ��Ҳ��������������Ϊ��Ӧֵ�������ԭ���ߵ�˼· ' , '
%           ʹ��ͨ������ֱֵ�ӻ�õģ���һ������Ӧֵ ' ��
%          ComputeFlagOfMutation = { 'by prob' , 'by steady' , 'by constant' } . �ֱ��ʾ�� '���������ֲ��������������ı仯�й�ʱ��������Ⱥ�ı�����' ,
%          '���������ֲ����ʽ���Ǩ��Ǩ�����йض����������ı仯�޹�ʱ��������Ⱥ�ı�����' ,
%          '���������ֲ����ʳʹ̶�ֵʱ��������Ⱥ�ı�����' ��
%         DisplayFlag = true or false, whether or not to display and plot results.
%         ProbFlag = true or false, whether or not to use probabilities to update emigration rates.
%         RandSeed = random number seed
% OUTPUTS: MinCost = array of best solution, one element for each generation
%          Hamming = final Hamming distance between solutions
%          �����֮���Hamming����
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

OPTIONS.pmodify = 1; % habitat modification probability ��Ϣ�ص�Ǩ�Ƹ���

OPTIONS.pmutate = 0.001; % initial mutation probability��ʼ��ͻ�����

Keep = 2; % elitism parameter: how many of the best habitats to keep from one generation to the next ��Ӣ����
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
%   ��ʼ��ÿһ����Ϣ�ص���Ⱥ�����ĸ���
for j = 1 : length(Population)
    Prob(j) = 1 / length(Population);
end




%% Begin the optimization loop
for GenIndex = 1 : OPTIONS.Maxgen   % �ӵ�һ����ʼ
    %% Save the best habitats in a temporary array.
    for j = 1 : Keep
        chromKeep(j,:) = Population(j).chrom;
        costKeep(j) = Population(j).cost;
    end
    
    %% ������Ӧֵ���� Map cost values to species counts.
    [Population] = GetFitnessCost(ComputeFlagOfFitnessFunction,Population, P);
    
    %% Migration ѡȡǨ����ģ�ͣ�����Ǩ����
    % Compute immigration rate and emigration rate for each species count.
    % lambda(i) is the immigration rate for habitat i.
    % mu(i) is the emigration rate for habitat i.
    [lambda, mu] = GetLambdaMu(ComputeFlagOfMigration,Population, I, E, P);
    
    
    %% �������ֲ�ͬ����Ⱥ�ֲ����ʼ��㷽ʽ
    % Compute the time derivative of Prob(i) for each habitat i.
    % ����ÿ����Ϣ�صĸ��ʶ�ʱ��ĵ���
    % ע�⣺�ٴ�ǿ�������к���ǰ���Ѿ��ٶ���Ϣ�� i�Ѿ����ŵ���������ˡ�
    switch ComputeFlagOfMutation         % ���ݲ����Ӧֵk����Ǩ���� If allowed to use probablities to update emigration rates .
        case 'by prob'         % ����ÿ����Ϣ�صĸ��ʶ�ʱ��ĵ���
            Prob=GetProbablisticProb(Population,Prob,I,E,P,lambda,mu,dt);
        case 'by steady'     % Ǩ������Ϊ����
            Prob=GetSteadyStateProb (Population,Prob,I,lambda,mu);
    end
    
    %% ���� BBO Ǩ�Ʒ�ʽ�����ӣ�
    [Island]=MigrationMethod(OPTIONS, OperatorFlagOfMigration, Population, lambda, mu, lambdaLower, lambdaUpper, muLower, muUpper);
    
    %%  Mutation ���������
    [Island]=ComputeMutation(OPTIONS, ComputeFlagOfMutation, Population, Island, Prob,MinParValue, MaxParValue);
    
    %% ����
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
    
    %% ��ʾ�����
    if DisplayFlag
        disp(['The best and mean of Generation # ', num2str(GenIndex), ' are ',...
            num2str(MinCost(end)), ' and ', num2str(AvgCost(end))]);
    end
    ExportData.BestAndMeanOfGeneration(GenIndex+1).generation(SET.monte_index)=GenIndex;
    ExportData.BestAndMeanOfGeneration(GenIndex+1).best(SET.monte_index)=MinCost(end);
    ExportData.BestAndMeanOfGeneration(GenIndex+1).mean(SET.monte_index)=AvgCost(end);

end

%% ����
[ExportData]=Conclude(SET, OPTIONS, DisplayFlag, ExportData, Population, nLegal, MinCost);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���ԸĽ�
%% Plot some results
% PlotConcludeFigures( SET, OPTIONS , 'BBO' , MinCost)

    
%% Obtain a measure of population diversity
% ����Ϣ��k����һ�㼶����ı���
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
%% GetFitnessCost ����
function [Population] = GetFitnessCost(ComputeFlagOfFitnessFunction,Population, P)
% Map cost values to fitness values.
% This loop assumes the population is already sorted from most fit to least fit.
% ѭ��ǰ�������ֵ�Ѿ��Ӻõ����ź���
switch ComputeFlagOfFitnessFunction
    % Map cost values to species counts.
    % ʹ�ò����Ӧֵ��Ҳ��������������Ϊ��Ӧֵ�������ԭ���ߵ�˼·
    case 'species count fitness'
        for i = 1 : length(Population)
            if Population(i).cost < inf
                Population(i).FitnessCost = P - i;
            else
                Population(i).FitnessCost = 0;
            end
        end
        
        % ʹ��ͨ������ֱֵ�ӻ�õģ���һ������Ӧֵ
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
%% GetLambdaMu����
function [lambda, mu] = GetLambdaMu(ComputeFlagOfMigration,Population, I, E, P)
switch ComputeFlagOfMigration
    case 'L-L'
        % Ǩ��ģ�Ͳ��� ����Ǩ��-����Ǩ�� ��ʽ
        % Compute immigration rate and extinction rate for each species count.
        % lambda(i) is the immigration rate for individual i.
        % mu(i) is the extinction rate for individual i.
        
        for i = 1 : length(Population)
            lambda(i) = I * (1 - Population(i).FitnessCost / P);
            mu(i) = E * Population(i).FitnessCost / P;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%% ���Ժ�����
        
end

return;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GetSteadyStateOfProb ����
function [Prob]=GetSteadyStateProb (Population,Prob,I,lambda,mu)
% ����Ǩ��Ǩ�����йض����������ı仯�޹�ʱ���������ֲ�����
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
%% GetProbablisticProb ����
function [Prob]=GetProbablisticProb(Population,Prob,I,E,P,lambda,mu,dt)
%% Compute the time derivative of Prob(i) for each habitat i.
% ����ÿ����Ϣ�صĸ��ʶ�ʱ��ĵ���
% ע�⣺�ٴ�ǿ�������к���ǰ���Ѿ��ٶ���Ϣ�� i�Ѿ����ŵ���������ˡ�
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
Prob = max(Prob, 0);    % ʹProb������ֵ�Ǹ�
Prob = Prob / sum(Prob);
return;
end


%% MigrationMethod Ǩ�����Ӻ���
function [Island]=MigrationMethod(OPTIONS, OperatorFlagOfMigration, Population, lambda, mu, lambdaLower, lambdaUpper, muLower, muUpper)
%% ���� BBO Ǩ�Ʒ�ʽ�����ӣ�
switch OperatorFlagOfMigration
    case 'PI'  % it means partial immigration-based BBO       
        %% Now use lambda and mu to decide how much information to share between habitats.
        % Ǩ�����㣺Ǩ���������ֱ���Ǩ�ơ��� lambda �� mu ȷ���� habitats ֮����Է������Ϣ����
        lambdaMin = min(lambda);
        lambdaMax = max(lambda);
        for k = 1 : length(Population)
            if rand > OPTIONS.pmodify
                continue;
            end
            % Normalize the immigration rate.
            % lamdaScale �����е� lambda(k) �ĳ߶�ѹ���� [0,1] ʵ���ռ䣬ʹ֮������rand()�Ƚ�
            lambdaScale = lambdaLower + (lambdaUpper - lambdaLower) * (lambda(k) - lambdaMin) ./ (lambdaMax - lambdaMin);
            %% Probabilistically input new information into habitat i
            for j = 1 : OPTIONS.numVar
                if rand < lambdaScale    % ����lambdaScaleȷ���������ֵ�Ƿ�Ǩ��
                    % Pick a habitat from which to obtain a feature
                    % ����mu ѡȡһ�� habitat ��ΪǨ���ߣ�ͨ�����̶�ѡȡ����muԽ���habitatԽ���ױ�ѡ��
                    RandomNum = rand * sum(mu);
                    Select = mu(1);
                    SelectIndex = 1;
                    while (Select < RandomNum) & (SelectIndex < OPTIONS.popsize)
                        SelectIndex = SelectIndex + 1;
                        Select = Select + mu(SelectIndex);
                    end
                    %����Ǩ�㣬����������ѡ�е���Ϣ�صĵ� j ��SIV���Ƹ���k����Ϣ��
                    Island(k,j) = Population(SelectIndex).chrom(j);
                else
                    Island(k,j) = Population(k).chrom(j);       % ������Ǩ��
                end
            end
        end
        

        
        
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %   �㷨����ʸ�����㷨��ʵ����������ٶȱ�����
%         %         ���룺
%         %             Population
%         %         ���̣�
%         %             �����ж�OPTIONS.pmodify��rand������randPmodify
%         %             ��OPTIONS.pmodify������������ʽ���������ж�lambda�Ƿ�Ǩ���rand������
%         %             Population����In
%         %             ��isPmodify��ÿ�������ж��Ƿ�Ǩ�ƣ�����Ǩ���߼�����L01����L01����In����Temp01��
%         %             ���ɱ�׼��Ǩ���ʾ���lambdaScale����=Ǩ����壬��=��������lambdaScaleȡ��lambda
%         %             ���������ж�lambda�Ƿ�Ǩ���rand������randMigration
%         %             ��lambdaScale�жϸ�����ֵ�Ƿ�Ǩ�ƣ�����Ǩ���߼�����L02������L02����In����Temp02��
%         %             ���̶�ѡ��һ��Ǩ�����壬����Index��
%         %                 �������̶�ָ��RandomPointer
%         %                 ������Ԫ�����mu���ۼƾ���Mu����͡����μ�ȥRandomPointer����������������1���������0��ת�����߼�����L03�����+1���õ�Ǩ�������±�Index��
%         %             Ǩ���������Ǩ�ƾ���Index��  ��=Ǩ����壬��=������
%         %             ��L01��L02����Index
%         %           ���������̶Ĺ��̲���Index����
%         %            ��Indexӳ�䵽Ǩ�ƹ�������Migration����=Ǩ����壬��=Ǩ�����壬ҳ=������
%         %           ��Migrationʵ��Ǩ������������Ǩ������Temp03��
%         %           ��Ǩ������Temp03��ֵ���� Island
%         %         �����
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
%         % ���̶�ѡ��һ��Ǩ������
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
%         % ����Index����
%         Migration=zeros(OPTIONS.popsize,OPTIONS.popsize,OPTIONS.numVar);
%         for i=1:OPTIONS.numVar
%             for j=1:OPTIONS.popsize
%                 if Index(j,i)==0
%                     continue;
%                 end
%                 Migration(j,Index(j,i),i)=1;
%             end
%         end
%         % ��ֵ����
%         for i=1:OPTIONS.numVar
%             Temp03(:,i)=Migration(:,:,i)*Temp02(:,i);
%         end
%             Island=In.chrom;
%             Island=Island.*~Temp03;
%             Island=Island+Temp03;
        

            
            
    case 'TI' % it means total immigration-based BBO
        %% Now use lambda and mu to decide how much information to share between habitats.
        % Ǩ�����㣺Ǩ������ȫ������Ǩ�ơ��� lambda �� mu ȷ���� habitats ֮����Է������Ϣ����
        lambdaMin = min(lambda);
        lambdaMax = max(lambda);
        for k = 1 : length(Population)
            if rand > OPTIONS.pmodify
                continue;
            end
            % Normalize the immigration rate.
            % lamdaScale �����е� lambda(k) �ĳ߶�ѹ���� [0,1] ʵ���ռ䣬ʹ֮������rand()�Ƚ�
            lambdaScale = lambdaLower + (lambdaUpper - lambdaLower) * (lambda(k) - lambdaMin) / (lambdaMax - lambdaMin);
            %% Probabilistically input new information into habitat i
            if rand < lambdaScale    % ����lambdaScaleȷ���������ֵ�Ƿ�Ǩ��
                for j=1 : OPTIONS.numVar   % ����ÿ������ֵ����Ҫѡ����Ե�Ǩ���߶�Ӧ������ֵ����Ǩ��
                    % ����mu ѡȡһ�� habitat ��ΪǨ���ߣ�ͨ�����̶�ѡȡ����muԽ���habitatԽ���ױ�ѡ��
                    RandomNum = rand * sum(mu);
                    Select = mu(1);
                    SelectIndex = 1;
                    while (Select < RandomNum) & (SelectIndex < OPTIONS.popsize)
                        SelectIndex = SelectIndex + 1;
                        Select = Select + mu(SelectIndex);
                    end
                    %% �������е� ����ֵ ������Ǩ��
                    %����Ǩ�㣬����������ѡ�е���Ϣ�صĵ� j ��SIV���Ƹ���k����Ϣ��
                    Island(k,j) = Population(SelectIndex).chrom(j);
                end
            else
                Island(k,:) = Population(k).chrom;       % ������Ǩ��
            end
        end
        
    case 'PE' % it means partial emigration-based BBO
        %% Now use lambda and mu to decide how much information to share between habitats.
        % Ǩ�����㣺Ǩ���������ֱ���Ǩ�ơ��� lambda �� mu ȷ���� habitats ֮����Է������Ϣ����
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
            % muScale �����е� lambda(k) �ĳ߶����ŵ� [0,1] ʵ���ռ䣬ʹ֮������rand()�Ƚ�
            muScale = muLower + (muUpper - muLower) * (mu(k) - muMin) / (muMax - muMin);
            %% Probabilistically input new information into habitat i
            for j = 1 : OPTIONS.numVar
                if rand < muScale    % ����lambdaScaleȷ���������ֵ�Ƿ�Ǩ��
                    % ���� lambda ѡȡһ�� habitat ��ΪǨ���ߣ�ͨ�����̶�ѡȡ����lambda Խ��� habitatԽ���ױ�ѡ��
                    RandomNum = rand * sum(lambda);
                    Select = lambda(1);
                    SelectIndex = 1;
                    while (Select < RandomNum) & (SelectIndex < OPTIONS.popsize)
                        SelectIndex = SelectIndex + 1;
                        Select = Select + lambda(SelectIndex);
                    end
                    %����Ǩ�㣬���������ѵ�k����Ϣ�صĵ� j ��SIV���Ƹ���ѡ�еĸ���Ϣ��
                    Island(SelectIndex,j) = Population(k).chrom(j);
                else
                    Island(k,j) = Population(k).chrom(j);       % ������Ǩ��
                end
            end
        end
        
    case 'TE' % it means total emigration-based BBO
        %% Now use lambda and mu to decide how much information to share between habitats.
        % Ǩ�����㣺Ǩ������ȫ������Ǩ�ơ��� lambda �� mu ȷ���� habitats ֮����Է������Ϣ����
        muMin=min(mu);
        muMax = max(mu);
        for ii=1:length(Population)
            Island(ii,:) = Population(ii).chrom;
        end
        for k=1:length(Population)
            if rand > OPTIONS.pmodify
                continue;
            end
            % ��muѡ��һ��������ΪǨ���ߣ�
            muScale = muLower+(muUpper-muLower)*(mu(k)-muMin)/(muMax-muMin);        % muScale �������и�������е� mu �߶ȵ� [lower,upper] �ռ䣬ʹ֮������rand()�Ƚ�
            if rand<muScale
                % �Խ�����ѡ�е���Ϣ�ص���������ֵ����Ǩ�Ʋ���
                for j=1 : OPTIONS.numVar
                    % �Ը�����ֵ����lambdaѡ��һ��������ΪǨ���ߡ������̶ĵķ�����lambda Խ��Խ���ױ�ѡ�С�
                    Select = lambda(1);
                    SelectIndex = 1;
                    RandomNum = rand * sum(lambda);
                    while (SelectIndex < OPTIONS.popsize ) & ( Select < RandomNum)
                        SelectIndex = SelectIndex+1;
                        Select = Select + lambda(SelectIndex);
                    end
                    % ����Ǩ�㣬�������ѵ�k����Ϣ�صĵ�j��SIV���Ƹ�ѡ�е���Ϣ��
                    Island(SelectIndex,j) = Population(k).chrom(j);
                end
            else
                Island(k,:) = Population(k).chrom;      % ������Ǩ��
            end
        end
        
end
return;
end

%%  ComputeMutation ��������ʺ���
function [Island]=ComputeMutation(OPTIONS, ComputeFlagOfMutation, Population, Island, Prob, MinParValue, MaxParValue)
switch ComputeFlagOfMutation
    case 'by prob'
        % ���������ֲ��������������ı仯�й�ʱ��������Ⱥ�ı�����
        Pmax = max(Prob);
        MutationRate = OPTIONS.pmutate * (1 - Prob / Pmax);
        
    case 'by steady'
        % ���������ֲ����ʽ���Ǩ��Ǩ�����йض����������ı仯�޹�ʱ��������Ⱥ�ı�����
        Pmax = max(Prob);
        MutationRate = OPTIONS.pmutate * (1 - Prob / Pmax);
        
    case 'by constant'
        % ���������ֲ����ʳʹ̶�ֵʱ��������Ⱥ�ı�����
        num=ones(1,OPTIONS.popsize)*OPTIONS.popsize;
        MutationRate=OPTIONS.pmutate * num;
        
end
% Mutate only the worst half of the solutions ���ñ���������һ�����Ⱥ����ͻ��
Population = PopSort(Population);
for index = round(length(Population)/2) : length(Population)    % ��ʱPopulation�Ѿ��ź�����
    for parnum = 1 : OPTIONS.numVar
        if rand < MutationRate(index)
            % �µ�ͻ��ֵ�����ɲ��ø�˹�ֲ������ɷ�ʽ
            Island(index,parnum) = floor(MinParValue + (MaxParValue - MinParValue + 1) * randn);
        end
    end
end

return;
end