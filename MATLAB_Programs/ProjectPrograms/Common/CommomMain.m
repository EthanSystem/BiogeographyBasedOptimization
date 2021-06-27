
%% Main 函数

clear;
clc;
diary off;

%% 建立数据文件夹
folderName01=([datestr(now,'yyyyMMddHHmmss')]);
folderParentPath01=(['D:\EthanLin\CoreFiles\ProjectsFile\Study\PostGraduate\Projects\BiogeographyBasedOptimization\MATLAB_Programs\ProjectPrograms\Common\Export Data']);
folderPath01=([folderParentPath01,'\',folderName01]);
sentence01=(['mkdir(''',folderPath01,''');']);
eval(sentence01);


sentence02=(['diary(''',folderPath01,'\CommandWindowRecord.txt'');']);
eval(sentence02);
diary on;


%% 计时设置


% tic;


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




%% OPTIONS 参数设置项
% WARNING: some of the optimization routines will not work if population size is odd.
MainOPTIONS.popsize = 8; % total population size 个体总数量
MainOPTIONS.Maxgen = 10; % generation count limit max number of gene  生成后代代数限制 
MainOPTIONS.numVar = 6; % number of genes in each population member   每个个体的基因数量
MainOPTIONS.pmutate = 0.01; % mutation probability     突变概率
MainOPTIONS.TimesOfMonte = 4; % times of Monte    蒙特卡罗模拟次数

%% 其它参数设置项
ExportData=0;
n=1;             % 重复次数，仅作为测试用，正常情况下=1



%% 选用的算法和测试基准函数
% Optimization methods
OptimizationFunction = [
% 'ACO   '; % ant colony optimization
'BBO   '; % biogeography-based optimization
% 'BBO_PI'; % biogeography-based optimization
% 'DE    '; % differential evolution
% 'ES    '; % evolutionary strategy
% 'GA    '; % genetic algorithm
% 'PBIL  '; % probability based incremental learning
% 'PSO   '; % particle swarm optimization
% 'StudGA'; % stud genetic algorithm
]

% Benchmark functions
BenchmarkFunction = [     %     multimodal? separable?  regular?   对应表格

'Sphere          '; %     n           y           y       1
% 'Step_ameliorate '; %     n           y           n       6
% 'Rosenbrock      '; %     n           n           y       5
% 'Quartic         '; %     n           y           y       7
% 'Ackley          '; %     y           n           y       10
% 'Griewank        '; %     y           n           y       11 
% 'Penalty1        '; %     y           n           y       12
% 'Penalty2        '; %     y           n           y       13
% 'Schwefel        '; %     y           y           n       8
% 'Rastrigin       '; %     y           y           y       9 
 ]



%% SET 参数设置

SET.row_of_subplot = size(OptimizationFunction,1);
SET.col_of_subplot = size(BenchmarkFunction,1);
SET.index_of_subplot = 1;
SET.monte_index=1;
SET.export_data_folder_path=folderPath01;


%% 主体部分：用来蒙特卡罗模拟，分别调用BBO的四种迁移算子

%% %%%测试程序运行速度的
for i=1:n
[ExportData] = Monte(SET ,MainOPTIONS, ExportData, OptimizationFunction , BenchmarkFunction , false);
end


%% 输出结果到 Export Data
sentence010=(['save(''',SET.export_data_folder_path,'\','Export Data'',''ExportData'');']);
eval(sentence010);

% 
% time01=toc;
% time01=time01/n;
% time01
% 


diary off;

%% 
sp=actxserver('SAPI.SpVoice');
text='程序运行完毕';
sp.Speak(text)
































%%