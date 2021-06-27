%% 这次测试的实验一 按照 文献 [4] 的实验二 进行

%% Copy of Optimization methods and Benchmark functions
% % Optimization methods
% OptFunction = [
% 'ACO   '; % ant colony optimization
% 'BBO   '; % biogeography-based optimization
% 'DE    '; % differential evolution
% 'ES    '; % evolutionary strategy
% 'GA    '; % genetic algorithm
% 'PBIL  '; % probability based incremental learning
% 'PSO   '; % particle swarm optimization
% 'StudGA']; % stud genetic algorithm
%
% % Benchmark functions
%  Bench = [     %     multimodal? separable?  regular?
%  'Ackley    '; %     y           n           y
%  'Fletcher  '; %     y           n           n
%  'Griewank  '; %     y           n           y
%  'Penalty1  '; %     y           n           y
%  'Penalty2  '; %     y           n           y
%  'Quartic   '; %     n           y           y
%  'Rastrigin '; %     y           y           y
%  'Rosenbrock'; %     n           n           y
%  'Schwefel  '; %     y           y           n
%  'Schwefel2 '; %     n           n           y
%  'Schwefel3 '; %     y           n           n
%  'Schwefel4 '; %     n           n           n
%  'Sphere    '; %     n           y           y
%  'Step      ']; %    n           y           n
%
% % for example : Bench = ['MAPSS'];


clear;
clc;





time01=tic

%%
TimesOfMonte = 4;






%%
% Optimization methods
OptimizationFunction = [
%     'BBO_PI   '; % partial immigration-based biogeography-based optimization
%     'BBO_TI   '; % total immigration-based biogeography-based optimization
    'BBO_PE   '; % partial emigration-based biogeography-based optimization
%     'BBO_TE   '; % total emigration-based biogeography-based optimization
    ];

% Benchmark functions
BenchmarkFunction = [     %     multimodal? separable?  regular?
    'Sphere    '; %     n           y           y
%     'Rosenbrock'; %     n           n           y
%     'Step      '  %    n           y           n
%     'Quartic   '; %     n           y           y
%     'Schwefel  '; %     y           y           n
%     'Rastrigin '; %     y           y           y
%     'Ackley    '; %     y           n           y
%     'Penalty1  '; %     y           n           y
%     'Penalty2  '; %     y           n           y
    ];



SET.row_of_subplot = size(OptimizationFunction,1);
SET.col_of_subplot = size(BenchmarkFunction,1);
SET.index_of_subplot = 1;



%% 主体部分：用来蒙特卡罗模拟，分别调用BBO的四种迁移算子
[MeanMin, MeanMinNorm, BestMin, BestMinNorm, MeanCPU] = Monte(SET , OptimizationFunction , BenchmarkFunction , TimesOfMonte , false)

time01=toc

sp=actxserver('SAPI.SpVoice');
text='程序运行完毕';
sp.Speak(text)
