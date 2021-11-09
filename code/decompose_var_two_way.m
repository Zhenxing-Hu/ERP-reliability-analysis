function [SStotal,var_trait, var_state, var_noise] = decompose_var_two_way(M)
    [n, k] = size(M);
    SStotal = var(M(:)) *(n*k - 1);
    MSR = var(mean(M, 2)) * k;
    MSW = sum(var(M,0, 2)) / n;
    MSC = var(mean(M, 1)) * n; 
    MSE = (SStotal - MSR *(n - 1) - MSC * (k -1))/ ((n - 1) * (k - 1));
    var_trait =  (MSR - MSE)/k;
    var_state = MSW - MSE;
    var_noise = MSE;
end