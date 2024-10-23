function x4 = computeeffectconc(BIS,Ce50,BIS_baseline,gamma)
% Compute effect site concentration for a given BIS (reverting hill
% function)
    x4 = Ce50.*(((BIS_baseline./BIS) -1) .^(1./gamma));
end