
%% Main ����
clear;
clc;
diary off;




%% �޸�����Ҫ�����
% numOfComputeFlagOfFitnessFunction=input('������BBO����Ӧֵ���㷽����1=species count fitnes��2=normalized fitness����');
% numOfComputeFlagOfMutation=input('������BBO�ı����ʼ��㷽����1=by prob��2=by steady��3=by constant����');

% for i=1:sizeOfExpressioin020(3)  ���¿���ɾ��
%     str030(i)=str010(1,:,numOfComputeFlag,numOfComputeFlagOfMutation);
% end




%% ���������ļ���
str020=datestr(now,'yyyymmdd');
folderName020=(str020);
str025='003';         % �޸�����Ҫ���õ����ļ��е�����
folderName025=str025;
% folderName030=str030;   ����ɾ��
folderParentPath01=(['D:\EthanLin\CoreFiles\ProjectsFile\Study\PostGraduate\Projects\BiogeographyBasedOptimization\MATLAB_Programs\ProjectPrograms\Project of BBO by HUNNU\Export Data']);
folderPath01=([folderParentPath01,'\',folderName020,'\',folderName025]);
sentence015=(['mkdir(''',folderPath01,''');']);
eval(sentence015);




%% OPTIONS ����������
% WARNING: some of the optimization routines will not work if population size is odd.
MainOPTIONS.popsize = 20; % total population size ����������
MainOPTIONS.Maxgen = 10; % generation count limit max number of gene  ���ɺ����������
MainOPTIONS.numVar = 6; % number of genes in each population member   ÿ������Ļ�������
MainOPTIONS.pmutate = 0.01; % mutation probability     ͻ�����
MainOPTIONS.TimesOfMonte = 3; % times of Monte    ���ؿ���ģ�����

%% ��������������
n=1;             % �ظ�����������Ϊ���Գ��������ٶ��ã������Ե������=1


%% ���� FunctionOptionsExpressions
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

%% �ֶ�����ѡ�õ��㷨�Ͳ��Ի�׼����
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
BenchmarkFunction = char({     %     multimodal? separable?  regular?   ��Ӧ���
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

%% ������Ҫ���㷨�Ĳ���ѡ�� FunctionOptionsExpressions
% ���� OptimizationFunction �� BBO �㷨�Ĳ���
for i=1:size(OptimizationFunction,1)
    string010(i,:)=strsplit(OptimizationFunction(i,:));
    string015(i,:)=char(string010(i,1));
    string020(i,:)=char(string010(i,end));
    EXPRESSIONS(i,:)=eval(string020(i,:));
end
OptimizationFunction=string015;     % ��������ȥ���˲�������Ż���������



% ���¿���ɾ��
% % �� expression010 �������ά�ϲ���һ��ά�ȣ��Ա����size��
% sizeOfExpression010=size(FunctionOptionsExpressions);
% expression020=reshape(FunctionOptionsExpressions,1,sizeOfExpression010(2),[]);
% sizeOfExpressioin020=size(expression020);








%% SET ��������

SET.row_of_subplot = size(OptimizationFunction,1);
SET.col_of_subplot = size(BenchmarkFunction,1);
SET.index_of_subplot = 1;
SET.monte_index=1;
SET.export_data_folder_path=folderPath01;
SET.numOfFunction=size(OptimizationFunction,1);      % ע�⣺�����������֮ǰ�� SET.numOfFunction=size(OptimizationFunction,1);
SET.numOfBench=size(BenchmarkFunction,1);
% SET.numOfFunctionOptions=sizeOfExpressioin020(3);    % ���������Ŀǰ�ķ�������ʱ�ò�����

%% д�� ExportData �Ļ�������
sentence020=(['ExportData.title='' '';']);
eval(sentence020);
% sentence020=(['ExportData.OPTIONS=',MainOPTIONS,';']);
% eval(sentence020);
% sentence030=(['ExportData.EXPRESSONS=',EXPRESSIONS,';']);
% eval(sentence030);


%% ��ʼ��¼�����������������
sentence050=(['diary(''',folderPath01,'\CommandWindowRecord.txt'');']);
eval(sentence050);
diary on;


%%%%%���Գ��������ٶȵ�ʱ����n����1��ֵ��
mainTimeStart=tic;
for i=1:n
    %% ���岿��
    % compute and display to window
    numOfFunction =SET.numOfFunction;
    numOfBench = SET.numOfBench;
%     numOfFunctionOptions = sizeOfExpressioin020(3);
    MeanMinCost = zeros(numOfFunction, numOfBench);
    BestMinCost = inf(numOfFunction, numOfBench);
    MeanCPUTime = zeros(numOfFunction, numOfBench);
    for indexFunction = 1 : numOfFunction       % �����Ż��㷨����
        SET.indexOfOptimizationFunction=indexFunction;
        for indexBench = 1 : numOfBench             % �������Ի�׼����
            SET.indexOfBenchmarkFunction=indexBench;
            disp(['Optimization method ', num2str(indexFunction), '/', num2str(numOfFunction), '  , Benchmark function ', num2str(indexBench), '/', num2str(numOfBench)]);
            %% �������ؿ���ģ��
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

mainTime=toc(mainTimeStart)/n;         % ���Գ��������ٶ��ã������ n ���ʹ��

%% ֹͣ��¼�����������������
diary off;

%% �������� Export Data
sentence060=(['save(''',SET.export_data_folder_path,'\','Export Data'',''ExportData'');']);
eval(sentence060);




%% �����������
sp=actxserver('SAPI.SpVoice');
text='�����������';
sp.Speak(text)






%% ����
% %% ѡ�õ��㷨�Ͳ��Ի�׼����
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
% BenchmarkFunction = [     %     multimodal? separable?  regular?   ��Ӧ���
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


