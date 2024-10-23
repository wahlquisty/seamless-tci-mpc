function r = get_refconc(N, E0, Ce50, gamma)
% Get reference concentration vector r for a given E0, Ce50 and gamma
    Ce_BIS50 = Ce50*((50-E0)/(E0-50-E0))^(1/gamma);
    r = Ce_BIS50 * ones(N,1);
end