
clear;
clc;


load('D:\EthanLin\CoreFiles\ProjectsFile\Study\PostGraduate\Projects\BiogeographyBasedOptimization\MATLAB_Programs\ProjectPrograms\Project of BBO by HUNNU\Export Data\20160519\017\Export Data.mat')


% 把曲线图画出来。
num_of_generation=length(ExportData.Details(1,1).BestAndMeanInEachGenerationAndMonte);
num_of_optFunction=size(ExportData.Details,1);
num_of_benchmarkFunction=size(ExportData.Details,2);
num_of_monte=length(ExportData.Details(1,1).RandomSeed);


for i=1:num_of_benchmarkFunction
    hold off;
    for f=1:num_of_optFunction
        for m=1:num_of_monte
            for g=1:num_of_generation
                Generation(g) = ExportData.Details(f,i).BestAndMeanInEachGenerationAndMonte(g).generation(1);
                eachBest(m,g) = ExportData.Details(f,i).BestAndMeanInEachGenerationAndMonte(g).best(m);
                eachMean(m,g) = ExportData.Details(f,i).BestAndMeanInEachGenerationAndMonte(g).mean(m);
            end
        end
        
        Best=mean(eachBest);
        Mean=mean(eachMean);
        
        % 画折线图
        subplot(num_of_benchmarkFunction,1,i);
        x=0:1:Generation(end);
        semilogy(x,Best);
        plot(x,Best,'-');
        hold on;
        
       % 画箱线图
       
    end
    hold off;
end











% %% 代码备份
% % 把曲线图画出来。
% num_of_generation=length(ExportData.Details(1,1).BestAndMeanOfGeneration);
% num_of_optFunction=size(ExportData.Details,1);
% num_of_benchmarkFunction=size(ExportData.Details,2);
% num_of_monte=length(ExportData.Details(1,1).RandomSeed);
% 
% 
% for i=1:num_of_benchmarkFunction
%     hold off;
%     for f=1:num_of_optFunction
%         for m=1:num_of_monte
%             for g=1:num_of_generation
%                 Generation(g) = ExportData.Details(f,i).BestAndMeanOfGeneration(g).generation(1);
%                 eachBest(m,g) = ExportData.Details(f,i).BestAndMeanOfGeneration(g).best(m);
%                 eachMean(m,g) = ExportData.Details(f,i).BestAndMeanOfGeneration(g).mean(m);
%             end
%         end
%         
%         Best=mean(eachBest);
%         Mean=mean(eachMean);
%         
%         % 画折线图
%         subplot(num_of_benchmarkFunction,1,i);
%         x=0:1:Generation(end);
%         semelogy()
%         plot([0:Generation(end)],Best,'-');
%         hold on;
%         
%        % 画箱线图
%        
%     end
%     hold off;
% end
% 
% 









