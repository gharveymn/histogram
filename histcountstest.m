function res = histcountstest(ts)
	
	for i = 1:numel(ts.tests)
		args = {};
		if(~isempty(ts.tests(i).BinEdges))
			args = [args, {'BinEdges', ts.tests(i).BinEdges}];
		end
		
		if(~isempty(ts.tests(i).NumBins))
			args = [args, {'NumBins', ts.tests(i).NumBins}];
		end
		
		if(~isempty(ts.tests(i).Normalization))
			args = [args, {'Normalization', ts.tests(i).Normalization}];
		end
		
		if(~isempty(ts.tests(i).BinMethod))
			args = [args, {'BinMethod', ts.tests(i).BinMethod}];
		end
		
		[res(i).c, res(i).e, res(i).b] = histcounts(ts.data, args{:});
		
	end
	
end