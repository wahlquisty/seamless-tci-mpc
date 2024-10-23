function BIS = computeBIS(x, BIS_baseline, Ce50, gamma)
% Calculate BIS from concentration
    BIS = BIS_baseline .* (Ce50 .^ gamma ./ (Ce50 .^ gamma + x .^ gamma));
end