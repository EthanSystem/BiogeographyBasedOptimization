function [ExportData]=Conclude(SET, OPTIONS, DisplayFlag, ExportData, Population, nLegal ,MinCost)
% Output results of population-based optimization algorithm.


if DisplayFlag
    % Count the number of duplicates
    NumDups = 0;
    for i = 1 : OPTIONS.popsize
        Chrom1 = sort(Population(i).chrom);
        for j = i+1 : OPTIONS.popsize
            Chrom2 = sort(Population(j).chrom);
            if isequal(Chrom1, Chrom2)
                NumDups = NumDups + 1;
            end
        end
    end
    disp([num2str(NumDups), ' duplicates in final population.']);
    ExportData.Details(SET.indexOfOptimizationFunction,SET.indexOfBenchmarkFunction).NumOfDuplicatesInFinalPopulation(SET.monte_index)=NumDups;
    
    disp([num2str(nLegal), ' legal individuals in final population.']);
    ExportData.Details(SET.indexOfOptimizationFunction,SET.indexOfBenchmarkFunction).NumOfLegalInFinalPopulation(SET.monte_index)=nLegal;

    % Display the best solution
    Chrom = sort(Population(1).chrom);
    disp(['Best chromosome = ', num2str(Chrom)]);
    ExportData.Details(SET.indexOfOptimizationFunction,SET.indexOfBenchmarkFunction).BestChromsome(:,SET.monte_index)=Chrom;

end


return;

end