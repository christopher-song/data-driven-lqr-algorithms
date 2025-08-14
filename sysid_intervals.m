%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% method for system ID when the number of intervals is varied
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rks_id, K_errs_id, P_errs_id] = sysid_intervals(intervals, dim_x, dim_u, x_diffs, xs, us, Q, R, K_true, P_true)
    % arrays for K and P and their errors
    K_errs_id = zeros(intervals,1);
    P_errs_id = zeros(intervals,1);
    rks_id = zeros(intervals,1);
    % do system ID with first i intervals of data
    for i = 1:intervals
        [rk, K_err, P_err] = sysid(dim_x, dim_u, x_diffs(1:i,:), xs(1:i,:), us(1:i,:), Q, R, K_true, P_true);
        % save data
        K_errs_id(i,1) = K_err;
        P_errs_id(i,1) = P_err;
        rks_id(i,1) = rk;
    end
end