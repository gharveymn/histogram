function res = histcountstest(ts)
	
	res = cell(numel(ts), 1);
	for j = 1:numel(ts)
		currts = ts(j);
		res{j}(numel(ts(j)), 4) = struct('c', [], 'e', [], 'b', []);
		
		for i = 1:numel(currts.tests)
			numargs = 0;
			args = cell(1,14);
			
			if(~isempty(currts.tests(i).NumBinsRaw))
				args{numargs+1} = currts.tests(i).NumBinsRaw;
				numargs = numargs+1;
			end
			
			if(~isempty(currts.tests(i).BinEdgesRaw))
				args{numargs+1} = currts.tests(i).BinEdgesRaw;
				numargs = numargs+1;
			end
			
			if(~isempty(currts.tests(i).NumBins))
				args{numargs+1} = 'NumBins';
				args{numargs+2} = currts.tests(i).NumBins;
				numargs = numargs+2;
			end
			
			if(~isempty(currts.tests(i).BinEdges))
				args{numargs+1} = 'BinEdges';
				args{numargs+2} = currts.tests(i).BinEdges;
				numargs = numargs+2;
			end

			if(~isempty(currts.tests(i).BinLimits))
				args{numargs+1} = 'BinLimits';
				args{numargs+2} = currts.tests(i).BinLimits;
				numargs = numargs+2;
			end
			
			if(~isempty(currts.tests(i).BinWidth))
				args{numargs+1} = 'BinWidth';
				args{numargs+2} = currts.tests(i).BinWidth;
				numargs = numargs+2;
			end
			
			if(~isempty(currts.tests(i).Normalization))
				args{numargs+1} = 'Normalization';
				args{numargs+2} = currts.tests(i).Normalization;
				numargs = numargs+2;
			end

			if(~isempty(currts.tests(i).BinMethod))
				args{numargs+1} = 'BinMethod';
				args{numargs+2} = currts.tests(i).BinMethod;
				numargs = numargs+2;
			end
			
			% do tests with 0, 1, 2, and 3 returns
			histcounts(currts.data, args{1:numargs});
			res{j}(i,1).c = ans;
			
			histcounts(currts.data, args{1:numargs});
			[res{j}(i,2).c] = histcounts(currts.data, args{1:numargs});
			
			histcounts(currts.data, args{1:numargs});
			[res{j}(i,3).c, res{j}(i,3).e] = histcounts(currts.data, args{1:numargs});
			
			[res{j}(i,4).c, res{j}(i,4).e, res{j}(i,4).b] = histcounts(currts.data, args{1:numargs});
		end
	end
end