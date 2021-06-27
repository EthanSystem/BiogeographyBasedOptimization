function [Population] = ClearDups(Population, MaxParValue, MinParValue)

% Make sure there are no duplicate individuals in the population.
% This logic does not make 100% sure that no duplicates exist, but any duplicates that are found are
% randomly mutated, so there should be a good chance that there are no duplicates after this procedure.
for i = 1 : length(Population)
	Chrom1 = sort(Population(i).chrom);
	for j = i+1 : length(Population)
		Chrom2 = sort(Population(j).chrom);
		if isequal(Chrom1, Chrom2)
            % 选取染色体上的某个位，对该位的值赋新值。这样做为了在这个函数的处理后，所有个体的染色体尽量不出现重复，当然难以保证可能还会出现重复值。
			parnum = ceil(length(Population(j).chrom) * rand);        
			Population(j).chrom(parnum  ) = floor(MinParValue + (MaxParValue - MinParValue + 1) * rand);
		end
	end
end
return;
