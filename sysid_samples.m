%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% method for system ID when the number of samples per interval is varied
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rks_id, K_errs_id, P_errs_id] = sysid_samples(nsamples, dim_x, dim_u, x_diffs_sampled, xs_sampled, us_sampled, Q, R, K_true, P_true)
    % arrays for K and P and their errors
    K_errs_id = zeros(nsamples,1);
    P_errs_id = zeros(nsamples,1);
    rks_id = zeros(nsamples,1);
    % do system ID with first i intervals of data
    for i = 1:nsamples
        [rk, K_err, P_err] = sysid(dim_x, dim_u, x_diffs_sampled(:,:,i), xs_sampled(:,:,i), us_sampled(:,:,i), Q, R, K_true, P_true);
        % save data
        K_errs_id(i,1) = K_err;
        P_errs_id(i,1) = P_err;
        rks_id(i,1) = rk;
    end
end