function entropy = getEntropy(Population, ub, lb, m)

    [dec, ~] = getPopulationDecAndObj(Population);
    X = sort(dec);
    N = size(X,1);
    Q1 = X(int32(N*1/4),:);
    Q2 = X(int32(N*2/4),:);
    Q3 = X(int32(N*3/4),:);
    
    
    inf = (Q3 - Q1)./(ub - lb);
    inf = inf + 0.000001;
    lgInf = log10(inf);

    mid = abs((Q2 - m)./(ub - lb));
    mid = mid + 0.000001;
    lgMid = log10(mid);
    
    entropy = - sum(inf.*lgInf,2) - sum(mid.*lgMid,2);
    
%     entropy = - sum(inf.*lgInf,2);

end

