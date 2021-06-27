
%% Main 函数

clear;
clc;


%% 计时设置





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




%% 参数设置项，初始化

OPTIONS.DisplayFlag=false;
OPTIONS.popsize=100;
OPTIONS.numVar=64;
OPTIONS.randNum=zeros(OPTIONS.popsize,OPTIONS.numVar);





%% 选用的测试基准函数


% Benchmark functions
BenchmarkFunction = [     %     multimodal? separable?  regular?   对应表格
    %  'Sphere          '; %     n           y           y       1
    %  'Schwefel3       '; %     y           n           n       2
    %  'Schwefel2       '; %     n           n           y       3
    %  'Schwefel4       '; %     n           n           n       4
    %  'Rosenbrock      '; %     n           n           y       5
    'Step            '; %     n           y           n       6
    'Step_ameliorate '; %     n           y           n       6
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
    ]









%%  %%% 增加随机数矩阵，行=个体数，列=变量维度
for i=1:OPTIONS.popsize
    OPTIONS.randNum(i,:)=rand(1,OPTIONS.numVar);
end



%% %%% 计算原始测试函数的值和要验证的测试函数的值

expression01=(['ComputeBenchmarks(@',BenchmarkFunction(1,:),',OPTIONS);']);
[population01,cost01]=eval(expression01);


expression02=(['ComputeBenchmarks(@',BenchmarkFunction(2,:),',OPTIONS);']);
[population02,cost02]=eval(expression02);


%% %%%计算二者差距

minus=cost02-cost01;

distance=sum(minus)







%% 程序结束播音
sp=actxserver('SAPI.SpVoice');
text='程序运行完毕';
sp.Speak(text)




















