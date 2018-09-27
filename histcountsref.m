function tests = histcountsref(x)
	normtypes = {'count','probability','countdensity','pdf','cumcount','cdf'};
	bms       = {'auto','scott','fd','integers','sturges','sqrt'};

	tests = cell(numel(normtypes)*numel(bms),1);


	for i = 1:numel(normtypes)
		nt = normtypes{i};
		for j = 1:numel(bms)
			bm = bms{j};

			s.Normalization = nt;
			s.BinMethod     = bm;
			[s.counts,s.edges,s.binidx] = histcounts(x, 'Normalization',...
											 nt, 'BinMethod', bm);

			tests{(i-1)*numel(normtypes) + j} = s;
		end
	end
end