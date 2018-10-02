function res = histcountstest(ts)
	
	res = cell(1, numel(ts));
	for j = 1:numel(ts)
		currts = ts(j);
		for i = 1:numel(currts.tests)
			args = {};
			
			if(~isempty(currts.tests(i).NumBins))
				args = [args, {'NumBins', currts.tests(i).NumBins}];
			end
			
			if(~isempty(currts.tests(i).BinEdges))
				args = [args, {'BinEdges', currts.tests(i).BinEdges}];
			end

			if(~isempty(currts.tests(i).BinLimits))
				args = [args, {'BinLimits', currts.tests(i).BinLimits}];
			end
			
			if(~isempty(currts.tests(i).BinWidth))
				args = [args, {'BinWidth', currts.tests(i).BinWidth}];
			end
			
			if(~isempty(currts.tests(i).Normalization))
				args = [args, {'Normalization', currts.tests(i).Normalization}];
			end

			if(~isempty(currts.tests(i).BinMethod))
				args = [args, {'BinMethod', currts.tests(i).BinMethod}];
			end

			[res{j}(i).c, res{j}(i).e, res{j}(i).b] = histcounts(currts.data, args{:});
		end
	end
	
end