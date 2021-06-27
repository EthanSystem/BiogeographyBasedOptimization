function [EachCPUTime, EachCost, ExportData] = Monte(SET , OPTIONS, EXPRESSIONS, ExportData, OptFunction,Bench,DisplayFlag)

% ���ؿ���ģ��
% Monte Carlo execution of population-based optimization software
% INPUT TimesOfMonte
% ����ֵ��������ȷ����Ĭ��ֵΪ�����ڲ�ȷ����ֵ
% OUTPUT MeanMin is the mean of the best solution found. It is a
% nFunction x nBench array, where nFunction is the number of optimization
% functions that are used, and nBench is the number of benchmarks that
% are optimized.
% OUTPUT MeanMinNorm is MeanMin normalized to a minimum of 1 for each benchmark.
% OUTPUT BestMin is the best solution found by each optimization function
% for each benchmark.
% OUTPUT BestMinNorm is BestMin normalized to a minimum of 1 for each benchmark.
% OUTPUT MeanCPU is the mean CPU time required for each optimization function
% normalized to 1.




% Ĭ����ʾ������
if ~exist('DisplayFlag','var')
    DisplayFlag = true;
end

nMonte = OPTIONS.TimesOfMonte;    % number of Monte Carlo runs




% compute and display to window
for indexMonte = 1 : nMonte
    monteTimeStart=tic;
    if DisplayFlag
        disp(['This is the #',num2str(indexMonte),' time']);
    end
    SET.monte_index=indexMonte;
    
    expression01=([OptFunction,'(SET,OPTIONS,ExportData,@',Bench,',',EXPRESSIONS,');']);
    [Cost, ~, ExportData] = eval(expression01);
    
    
    % Calculate cost in each monte
    monteTime=toc(monteTimeStart);
    EachCPUTime(indexMonte)=monteTime;
    EachCost(indexMonte)=Cost(end);
    
end













