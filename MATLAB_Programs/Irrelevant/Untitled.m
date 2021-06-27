clear;
clc;

n=1000;

A=rand(n);


time=zeros(1:10);

% tic;
% for i=1:n
%     C=ones(n);
%     for j=1:n
%
%         C(i,j)=A(i,j);
%     end
% end
% time(1)=toc;
%
%
% tic;
% for i=1:n
%     for j=1:n
%         B(i,j)=A(i,j);
%     end
% end
% time(2)=toc;
%
% A=rand(n);
% C=ones(n);
%
% tic;
% C=A;
% time(3)=toc;




% %%
% tic;
% 
% for trequency=1:100
%     MinParValue=1;
%     Granularity = 0.1;
%     popsize = 50;
%     numVar=20;
%     MaxParValue = floor(1 + 2 * 5.12 / Granularity);
%     % Initialize population
%     for popindex = 1 : popsize
%         chrom = floor(MinParValue + (MaxParValue - MinParValue + 1) * rand(1,numVar));
%         Population(popindex).chrom = chrom;
%     end
%     for popindex = 1 : popsize
%         Population(popindex).cost = 0;
%         for i = 1 : numVar
%             gene = Population(popindex).chrom(i);
%             x = (gene - MinParValue) / (MaxParValue - MinParValue) * 2 * 100 - 100;
%             Population(popindex).cost = Population(popindex).cost + (floor(x+0.5))^2;
%         end
%     end
%     for i = 1 : popsize
%         for k = 1 :numVar
%             Population(i).chrom(k) = max(Population(i).chrom(k), MinParValue);
%             Population(i).chrom(k) = min(Population(i).chrom(k), MaxParValue);
%         end
%     end
% end
% 
% time(4)=toc;




%%
tic;

for trequency=1:100
    Granularity = 0.1;
    popsize = 50;
    numVar=20;
    MinParValue=1;
    MaxParValue = floor(1 + 2 * 5.12 / Granularity);
    % Initialize population
    gene= floor(MinParValue + (MaxParValue - MinParValue + 1) * rand(popsize,numVar));
    % Compute the cost of each member in Population
    x = (gene - MinParValue) / (MaxParValue - MinParValue) * 2 * 100 - 100;
    t=floor(x+0.5).^2;
    F=sum(t,2);
    for i=1:popsize
        Population(i).chrom=gene(i,:);
        Population(i).cost=F(i,:);
    end
    for i = 1 : popsize
            Population(i).chrom = max(Population(i).chrom, MinParValue);
            Population(i).chrom = min(Population(i).chrom, MaxParValue);      
    end
end



time(5)=toc;

%%

for i=1:length(time)
    %     disp(['time(',i,')= ',time(i)])
    time(i)
end