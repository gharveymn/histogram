function teststruct = histcountstestgen
	
	x = rand(100000,1);
	
	edges     = {[], sort(rand(10, 1)), sort(rand(1, 10))};
	nbins     = {[], 1, 100};
	normtypes = {'count','probability','countdensity','pdf','cumcount','cdf'};
	bms       = {'auto','scott','fd','integers','sturges','sqrt'};

	A1 = allcomb(edges, {[]}, normtypes, {[]}, {[]}, {[]}, {[]});
	A2 = allcomb({[]}, nbins, normtypes, bms, {[]}, {[]}, {[]});
	
	A = [A1; A2];
	
	t = cell2struct(A,{'BinEdges','NumBins','Normalization',...
		                  'BinMethod', 'c', 'e', 'b'}, 2);
				   
	teststruct.data  = x;
	teststruct.tests = t;
				   
end