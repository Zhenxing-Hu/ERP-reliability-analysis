function [SST,var_trait, var_state, var_noise] = decompose_var_two_way(M)
    [n, k] = size(M);
    SST = var(M(:)) *(n*k - 1);
    MSBS = var(mean(M, 2)) * k;
    MSBM =  var(mean(M, 1)) * n; 
    MSE = (SST - MSBS *(n - 1) - MSBM * (k -1))/ ((n - 1) * (k - 1));
    var_trait =  (MSBS - MSE)/k;
    var_state = (MSBM-MSE)/n;
    var_noise = MSE;
end
