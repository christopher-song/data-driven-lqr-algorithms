%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% method for adp when the number of intervals is varied
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rks_adp, K_errs_adp, P_errs_adp, max_iter_err] = adp_intervals(intervals, iterations, dim_x, dim_u, K0, kron_x_diffs, xxs, xus, Q, R, K_true, P_true)
    % arrays for K and P iterates and their errors
    K_iterrs_adp = zeros(intervals,1);
    P_iterrs_adp = zeros(intervals,1);
    K_errs_adp = zeros(intervals,1);
    P_errs_adp = zeros(intervals,1);
    rks_adp = zeros(intervals,1);
    
    % do adp iteration with first i intervals of data
    for i = 1:intervals
        [rk, K_iterr, P_iterr, K_err, P_err] = adp(dim_x, dim_u, K0, iterations, kron_x_diffs(1:i,:), xxs(1:i,:), xus(1:i,:), Q, R, K_true, P_true);
        % save data
        K_iterrs_adp(i,1) = K_iterr;
        P_iterrs_adp(i,1) = P_iterr;
        K_errs_adp(i,1) = K_err;
        P_errs_adp(i,1) = P_err;
        rks_adp(i,1) = rk;
    end
    
    max_iter_err = max(abs(K_iterrs_adp)); % show error