
%% Main ����
clear;
clc;
diary off;


%% ����Ҫ���е����
% expression11=(['''species count fitness'',''PI'',''L-L'',''by prob'',true']);
% expression12=(['''species count fitness'',''PI'',''L-L'',''by steady'',true']);
% expression13=(['''species count fitness'',''PI'',''L-L'',''by constant'',true']);
% expression21=(['''normalized fitness'',''PI'',''L-L'',''by prob'',true']);
% expression22=(['''normalized fitness'',''PI'',''L-L'',''by steady'',true']);
% expression23=(['''normalized fitness'',''PI'',''L-L'',''by constant'',true']);
% % �޸�����Ҫ�����
% EXPRESSIONS=expression13;
% 
%%%%% ���·����ᱨ����ʱ���á�
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
        expression010(1,:,i,j)=(['''',BBO.ComputeFlagOfFitnessFunction(i,:),''',''PI'',''L-L'',''',BBO.ComputeFlagOfMutation(j,:),''',true']);
        str010(1,:,i,j)=([BBO.ComputeFlagOfFitnessFunction(i,:),' and ',BBO.ComputeFlagOfMutation(j,:)]);
    end
end



% �޸�����Ҫ�����
% numOfComputeFlagOfFitnessFunction=input('������BBO����Ӧֵ���㷽����1=species count fitnes��2=normalized fitness����');
% numOfComputeFlagOfMutation=input('������BBO�ı����ʼ��㷽����1=by prob��2=by steady��3=by constant����');
numOfComputeFlagOfFitnessFunction=1;
numOfComputeFlagOfMutation=3;


EXPRESSIONS=expression010(1,:,numOfComputeFlagOfFitnessFunction,numOfComputeFlagOfMutation);
str030=str010(1,:,numOfComputeFlagOfFitnessFunction,numOfComputeFlagOfMutation);








%% ���������ļ���
str020=datestr(now,'yyyymmdd');
folderName020=(str020);
str025='001';         % �޸�����Ҫ���õ����ļ��е�����
folderName025=str025;
folderName030=str030;
folderParentPath01=(['D:\EthanLin\CoreFiles\ProjectsFile\Study\PostGraduate\Projects\BiogeographyBasedOptimization\MATLAB_Programs\ProjectPrograms\Project of BBO by HUNNU\Export Data']);
folderPath01=([folderParentPath01,'\',folderName020,'\',folderName025,'\',folderName030]);
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



%% ѡ�õ��㷨�Ͳ��Ի�׼����
% Optimization methods
OptimizationFunction = char({
    % 'ACO   '; % ant colony optimization
    'BBO'; % biogeography-based optimization
    % 'BBO_PI'; % biogeography-based optimization
    % 'DE    '; % differential evolution
    % 'ES    '; % evolutionary strategy
    % 'GA    '; % genetic algorithm
    % 'PBIL  '; % probability based incremental learning
    % 'PSO   '; % particle swarm optimization
    % 'StudGA'; % stud genetic algorithm
    });

% Benchmark functions
BenchmarkFunction = char({     %     multimodal? separable?  regular?   ��Ӧ���
    'Sphere',                         %     n           y           y       1
    'Step',                         %     n           y           n       6
    'Rosenbrock',              %     n           n           y       5
    'Quartic',                    %     n           y           y       7
    'Ackley',                     %     y           n           y       10
    'Griewank',                %     y           n           y       11
    'Penalty1',               %     y           n           y       12
    'Penalty2',                %     y           n           y       13
    'Schwefel',               %     y           y           n       8
    'Rastrigin',                  %     y           y           y       9
    });


%% SET ��������

SET.row_of_subplot = size(OptimizationFunction,1);
SET.col_of_subplot = size(BenchmarkFunction,1);
SET.index_of_subplot = 1;
SET.monte_index=1;
SET.export_data_folder_path=folderPath01;



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

%% ���岿�֣��������ؿ���ģ��
%%%%%���Գ��������ٶȵ�ʱ����n����1��ֵ��

mainTimeStart=tic;

for i=1:n
    [MeanMinCost, BestMinCost, MeanCPUTime, MeanMinCostNorm, BestMinCostNorm, MeanCPUTimeNorm, ExportData] = Monte(SET ,MainOPTIONS, EXPRESSIONS, ExportData, OptimizationFunction , BenchmarkFunction , false);
end

mainTime=toc(mainTimeStart)/n;         % ���Գ��������ٶ��ã������ n ���ʹ��

%% ֹͣ��¼�����������������
diary off;

%% �������� Export Data
sentence060=(['save(''',SET.export_data_folder_path,'\','Export Data ','kkk',''',''ExportData'');']);
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


