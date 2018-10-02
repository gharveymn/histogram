function ts = histcountstestgen
	
	%% 1
	
	x = rand(100000,1);
	
	nbins     = {[], 1, 100};
	edges     = {[], sort(rand(10, 1)), sort(rand(1, 10)),...
			   sort([rand(50,1); 0.4213; 0.4213]),...
			   sort([rand(50,1); 0.4213; 0.4213; 0.4213])};
	binlims   = {[], [-1,2], [.3;.7]};
	binwidth  = {[], 5, 1, .3};
	normtypes = {[], 'count','probability','countdensity','pdf','cumcount','cdf'};
	bms       = {[], 'auto','scott','fd','integers','sturges','sqrt'};

	A1 = allcomb({[]}, edges, {[]}, {[]}, normtypes, {[]});
	A2 = allcomb(nbins, {[]}, binlims, binwidth, normtypes, bms);
	A = [A1; A2];
	
	t = cell2struct(A,{'NumBins','BinEdges','BinLimits', 'BinWidth',...
				    'Normalization','BinMethod'}, 2);
				   
	ts(1).data  = x;
	ts(1).tests = t;
	
	%% 2
	
	x = randi(50,[1000, 1]);
	
	nbins     = {[], 1, 100};
	edges     = {[], sort(randi(50, [70,1])), sort(randi(50, [1,70])),...
			   sort([randi(23,[16,1]); 5; 5]),...
			   sort([randi(23,[16,1]); 7; 7; 7]),...
			   sort(randi(500,[30,1]))};
	binlims   = {[], [-50,100], [13;20]};
	binwidth  = {[], 100, 7, 1};
	
	A1 = allcomb({[]}, edges, {[]}, {[]}, normtypes, {[]});
	A2 = allcomb(nbins, {[]}, binlims, binwidth, normtypes, bms);
	A = [A1; A2];
	
	t = cell2struct(A,{'NumBins','BinEdges','BinLimits', 'BinWidth',...
				    'Normalization','BinMethod'}, 2);
				   
	ts(2).data = x;
	ts(2).tests = t;
	
	%% 3
	
	bms       = {[], 'auto','scott','integers','sturges','sqrt'};
	
	A1 = allcomb({[]}, edges, {[]}, {[]}, normtypes, {[]});
	A2 = allcomb(nbins, {[]}, binlims, binwidth, normtypes, bms);
	A = [A1; A2];
	
	t = cell2struct(A,{'NumBins','BinEdges','BinLimits', 'BinWidth',...
				    'Normalization','BinMethod'}, 2);
	
	ts(3).data = uint8(x);
	ts(3).tests = t;
	
	%% 4
	
	nbins     = {[], 1, 100};
	edges     = {[], uint8(sort(randi(50, [70,1]))), uint8(sort(randi(50, [1,70]))),...
			   uint8(sort([randi(23,[16,1]); 5; 5])),...
			   uint8(sort([randi(23,[16,1]); 7; 7; 7])),...
			   uint8(sort(randi(500,[30,1])))};
	binlims   = {[], int8([-50,100]), uint8([13;20])};
	binwidth  = {[], uint8(100), uint8(7), uint8(1)};
	
	bms       = {[], 'auto','scott','fd','integers','sturges','sqrt'};
	
	A1 = allcomb({[]}, edges, {[]}, {[]}, normtypes, {[]});
	A2 = allcomb(nbins, {[]}, binlims, binwidth, normtypes, bms);
	A = [A1; A2];
	
	t = cell2struct(A,{'NumBins','BinEdges','BinLimits', 'BinWidth',...
				    'Normalization','BinMethod'}, 2);
				   
	ts(4).data = x;
	ts(4).tests = t;
	
	%% 5
	
	bms       = {[], 'auto','scott','integers','sturges','sqrt'};
	
	A1 = allcomb({[]}, edges, {[]}, {[]}, normtypes, {[]});
	A2 = allcomb(nbins, {[]}, binlims, binwidth, normtypes, bms);
	A = [A1; A2];
	
	t = cell2struct(A,{'NumBins','BinEdges','BinLimits', 'BinWidth',...
				    'Normalization','BinMethod'}, 2);
	
	ts(5).data = uint8(x);
	ts(5).tests = t;
	
	%% 6
	
	x = rand(300);
	
	nbins     = {[], 1, 100};
	edges     = {[], sort(rand(10, 1)), sort(rand(1, 10)),...
			   sort([rand(50,1); 0.4213; 0.4213]),...
			   sort([rand(50,1); 0.4213; 0.4213; 0.4213])};
	binlims   = {[], [-1,2], [.3;.7]};
	binwidth  = {[], 5, 1, .3};
	normtypes = {[], 'count','probability','countdensity','pdf','cumcount','cdf'};
	bms       = {[], 'auto','scott','fd','integers','sturges','sqrt'};

	A1 = allcomb({[]}, edges, {[]}, {[]}, normtypes, {[]});
	A2 = allcomb(nbins, {[]}, binlims, binwidth, normtypes, bms);
	A = [A1; A2];
	
	t = cell2struct(A,{'NumBins','BinEdges','BinLimits', 'BinWidth',...
				    'Normalization','BinMethod'}, 2);
				   
	ts(1).data  = x;
	ts(1).tests = t;
	
	
	
end