
%% Main 函数
clear;
clc;
diary off;




%% 修改你需要的语句
% numOfComputeFlagOfFitnessFunction=input('请输入BBO的适应值计算方法（1=species count fitnes，2=normalized fitness）：');
% numOfComputeFlagOfMutation=input('请输入BBO的变异率计算方法（1=by prob，2=by steady，3=by constant）：');

% for i=1:sizeOfExpressioin020(3)  以下可以删除
%     str030(i)=str010(1,:,numOfComputeFlag,numOfComputeFlagOfMutation);
% end




%% 建立数据文件夹
str020=datestr(now,'yyyymmdd');
folderName020=(str020);
str025='003';         % 修改你需要放置的子文件夹的名称
folderName025=str025;
% folderName030=str030;   可以删除
folderParentPath01=(['D:\EthanLin\CoreFiles\ProjectsFile\Study\PostGraduate\Projects\BiogeographyBasedOptimization\MATLAB_Programs\ProjectPrograms\Project of BBO by HUNNU\Export Data']);
folderPath01=([folderParentPath01,'\',folderName020,'\',folderName025]);
sentence015=(['mkdir(''',folderPath01,''');']);
eval(sentence015);




%% OPTIONS 参数设置项
% WARNING: some of the optimization routines will not work if population size is odd.
MainOPTIONS.popsize = 20; % total population size 个体总数量
MainOPTIONS.Maxgen = 10; % generation count limit max number of gene  生成后代代数限制
MainOPTIONS.numVar = 6; % number of genes in each population member   每个个体的基因数量
MainOPTIONS.pmutate = 0.01; % mutation probability     突变概率
MainOPTIONS.TimesOfMonte = 3; % times of Monte    蒙特卡罗模拟次数

%% 其它参数设置项
n=1;             % 重复次数，仅作为测试程序运行速度用，不测试的情况下=1


%% 生成 FunctionOptionsExpressions
BBO.ComputeFlagOfFitnessFunction=char({
    'species count fitness',
    'normalized fitness',
    });

BBO.OperatorFlagOfMigration=char({
    'PI',
    'TI',
    'PE',
    'TE',
    });
BBO.ComputeFlagOfMigration=char({
    'L-L',
    });
BBO.ComputeFlagOfMutation=char({
    'by prob',
    'by steady',
    'by constant',
    });


for i=1:size(BBO.ComputeFlagOfFitnessFunction,1)
    for j=1:size(BBO.ComputeFlagOfMutation,1)
        FunctionOptionsExpressions(1,:,i,j)=(['''',BBO.ComputeFlagOfFitnessFunction(i,:),''',''PI'',''L-L'',''',BBO.ComputeFlagOfMutation(j,:),''',true']);
        str010(1,:,i,j)=([BBO.ComputeFlagOfFitnessFunction(i,:),' and ',BBO.ComputeFlagOfMutation(j,:)]);
    end
end

%% 手动更改选用的算法和测试基准函数
% Optimization methods
OptimizationFunction = char({
    %     'ACO', % ant colony optimization
%     'BBO FunctionOptionsExpressions(1,:,1,1)',% biogeography-based optimization
%     'BBO FunctionOptionsExpressions(1,:,1,2)',% biogeography-based optimization
%     'BBO FunctionOptionsExpressions(1,:,1,3)',% biogeography-based optimization
%     'BBO FunctionOptionsExpressions(1,:,2,1)',% biogeography-based optimization
%     'BBO FunctionOptionsExpressions(1,:,2,2)',% biogeography-based optimization
    'BBO FunctionOptionsExpressions(1,:,2,3)',% biogeography-based optimization
    %     'DE', % differential evolution
    %     'ES', % evolutionary strategy
    %     'GA', % genetic algorithm
    %     'PBIL', % probability based incremental learning
    %     'PSO', % particle swarm optimization
    %     'StudGA', % stud genetic algorithm
    });

% Benchmark functions
BenchmarkFunction = char({     %     multimodal? separable?  regular?   对应表格
    'Sphere',                         %     n           y           y       1
%     'Step',                         %     n           y           n       6
%     'Rosenbrock',              %     n           n           y       5
%     'Quartic',                    %     n           y           y       7
    %     'Ackley',                     %     y           n           y       10
    %     'Griewank',                %     y           n           y       11
    %     'Penalty1',               %     y           n           y       12
    %     'Penalty2',                %     y           n           y       13
    %     'Schwefel',               %     y           y           n       8
    %     'Rastrigin',                  %     y           y           y       9
    });

%% 设置需要的算法的参数选项 FunctionOptionsExpressions
% 生成 OptimizationFunction 的 BBO 算法的参数
for i=1:size(OptimizationFunction,1)
    string010(i,:)=strsplit(OptimizationFunction(i,:));
    string015(i,:)=char(string010(i,1));
    string020(i,:)=char(string010(i,end));
    EXPRESSIONS(i,:)=eval(string020(i,:));
end
OptimizationFunction=string015;     % 重新生成去除了参数项的优化函数集合



% 以下可以删除
% % 把 expression010 的最后两维合并成一个维度，以便求出size。
% sizeOfExpression010=size(FunctionOptionsExpressions);
% expression020=reshape(FunctionOptionsExpressions,1,sizeOfExpression010(2),[]);
% sizeOfExpressioin020=size(expression020);








%% SET 参数设置

SET.row_of_subplot = size(OptimizationFunction,1);
SET.col_of_subplot = size(BenchmarkFunction,1);
SET.index_of_subplot = 1;
SET.monte_index=1;
SET.export_data_folder_path=folderPath01;
SET.numOfFunction=size(OptimizationFunction,1);      % 注意：玄机在这里。改之前是 SET.numOfFunction=size(OptimizationFunction,1);
SET.numOfBench=size(BenchmarkFunction,1);
% SET.numOfFunctionOptions=sizeOfExpressioin020(3);    % 这个，按照目前的方法，暂时用不到。

%% 写入 ExportData 的基本属性
sentence020=(['ExportData.title='' '';']);
eval(sentence020);
% sentence020=(['ExportData.OPTIONS=',MainOPTIONS,';']);
% eval(sentence020);
% sentence030=(['ExportData.EXPRESSONS=',EXPRESSIONS,';']);
% eval(sentence030);


%% 开始记录公共窗口输出的数据
sentence050=(['diary(''',folderPath01,'\CommandWindowRecord.txt'');']);
eval(sentence050);
diary on;


%%%%%测试程序运行速度的时候用n大于1的值。
mainTimeStart=tic;
for i=1:n
    %% 主体部分
    % compute and display to window
    numOfFunction =SET.numOfFunction;
    numOfBench = SET.numOfBench;
%     numOfFunctionOptions = sizeOfExpressioin020(3);
    MeanMinCost = zeros(numOfFunction, numOfBench);
    BestMinCost = inf(numOfFunction, numOfBench);
    MeanCPUTime = zeros(numOfFunction, numOfBench);
    for indexFunction = 1 : numOfFunction       % 遍历优化算法函数
        SET.indexOfOptimizationFunction=indexFunction;
        for indexBench = 1 : numOfBench             % 遍历测试基准函数
            SET.indexOfBenchmarkFunction=indexBench;
            disp(['Optimization method ', num2str(indexFunction), '/', num2str(numOfFunction), '  , Benchmark function ', num2str(indexBench), '/', num2str(numOfBench)]);
            %% 用来蒙特卡罗模拟
            [EachCPUTime, EachCost, ExportData] = Monte(SET, MainOPTIONS, EXPRESSIONS(indexFunction,:), ExportData, OptimizationFunction(indexFunction,:), BenchmarkFunction(indexBench,:), false);
            
            % subplot index
            SET.index_of_subplot = SET.index_of_subplot +1;
            % Calculate the mean CPU time , mean min cost , best min cost .
            MeanCPUTime(indexFunction,indexBench)=mean(EachCPUTime);
            MeanMinCost(indexFunction,indexBench)=mean(EachCost);
            BestMinCost(indexFunction,indexBench)=min(EachCost);
            % Export
            ExportData.Details(SET.indexOfOptimizationFunction,SET.indexOfBenchmarkFunction).MeanCPUTime(indexFunction,indexBench)=MeanCPUTime(indexFunction,indexBench);
            ExportData.Details(SET.indexOfOptimizationFunction,SET.indexOfBenchmarkFunction).MeanMinCost(indexFunction,indexBench)=MeanMinCost(indexFunction,indexBench);
            ExportData.Details(SET.indexOfOptimizationFunction,SET.indexOfBenchmarkFunction).BestMinCost(indexFunction,indexBench)=BestMinCost(indexFunction,indexBench);
        end
    end
    
    
    % Normalize the results
    if min(MeanMinCost) == 0
        MeanMinCostNorm = [];
    else
        MeanMinCostNorm = MeanMinCost * diag(1./min(MeanMinCost));
    end
    
    if min(BestMinCost) == 0
        BestMinCostNorm = [];
    else
        BestMinCostNorm = BestMinCost * diag(1./min(BestMinCost));
    end
    
    MeanCPUTime = min(MeanCPUTime);
    MeanCPUTimeNorm = MeanCPUTime * diag(1./min(MeanCPUTime));
    
    % Export
    ExportData.MeanMinCostNorm = MeanMinCostNorm;
    ExportData.BestMinCostNorm = BestMinCostNorm;
    ExportData.MeanCPUTimeNorm = MeanCPUTimeNorm;
    
    
    
    
end

mainTime=toc(mainTimeStart)/n;         % 测试程序运行速度用，与参数 n 配合使用

%% 停止记录公共窗口输出的数据
diary off;

%% 输出结果到 Export Data
sentence060=(['save(''',SET.export_data_folder_path,'\','Export Data'',''ExportData'');']);
eval(sentence060);




%% 程序结束提醒
sp=actxserver('SAPI.SpVoice');
text='程序运行完毕';
sp.Speak(text)






%% 备用
% %% 选用的算法和测试基准函数
% % Optimization methods
% OptimizationFunction = [
% 'ACO   '; % ant colony optimization
% 'BBO   '; % biogeography-based optimization
% 'DE    '; % differential evolution
% 'ES    '; % evolutionary strategy
% 'GA    '; % genetic algorithm
% 'PBIL  '; % probability based incremental learning
% 'PSO   '; % particle swarm optimization
% 'StudGA'; % stud genetic algorithm
% ]
%
% % Benchmark functions
% BenchmarkFunction = [     %     multimodal? separable?  regular?   对应表格
%  'Sphere          '; %     n           y           y       1
%  'Schwefel3       '; %     y           n           n       2
%  'Schwefel2       '; %     n           n           y       3
%  'Schwefel4       '; %     n           n           n       4
%  'Rosenbrock      '; %     n           n           y       5
%  'Step            '; %     n           y           n       6
%  'Quartic         '; %     n           y           y       7
%  'Schwefel        '; %     y           y           n       8
%  'Rastrigin       '; %     y           y           y       9
%  'Ackley          '; %     y           n           y       10
%  'Griewank        '; %     y           n           y       11
%  'Penalty1        '; %     y           n           y       12
%  'Penalty2        '; %     y           n           y       13
%  'Shekel          '; %     y                               14
%  'Kowalik         '; %     y                               15
%  'SixHumpCamelBack'; %     y                               16
%  'BraninRCOS      '; %     y                               17
%  'GoldsteinPrice  '; %     y                               18
%  'HartmanFamily19 '; %     y                               19
%  'HartmanFamily20 '; %                                     20
%  'Shekel21        '; %     y                               21
%  'Shekel22        '; %     y                               22
%  'Shekel23        '; %     y                               23
%  'Schaffer        '; %     y                               24
%  'Bohachevsky02   '; %                                     25
%  'Bohachevsky     '; %                                     26
%  'Schwefel        '; %                                     27
%  'Fletcher        ';%     y           n           n       28
%  ]


