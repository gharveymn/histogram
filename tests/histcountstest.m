function res = histcountstest(ts)
	res = cell(numel(ts), 1);
	for j = 1:numel(ts)
		currts = ts{j};
		res{j}(numel(currts.tests), 4) = struct('c', [], 'e', [], 'b', []);
		for i = 1:numel(currts.tests)
			res{j}(i,:) = testwithopts (currts.data, currts.tests(i));
		end
	end
end