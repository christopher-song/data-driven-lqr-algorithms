clear all
rng(6,"twister");
warning('off')
format long 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% choose system 2 or 3 or 6 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

system = 6;


nt = [3, 4, 5, 6, 7, 8, 10, 13, 16, 20, 25, 32, 40, 50, 63, 79, 100, 126, 158, 200, 251, 316, 398, 501, 631, 794, 1000, 1258, 1584, 1995, 2511, 3162, 3981, 5011, 6309, 7943, 10000]; % number of samples per interval
%nt = [10, 100, 1000];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters for simulation and data collection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if system == 2
    [Q, R, K0, A, B, a, f, dim_x, dim_u, K_true, P_true] = params_2d();
    x0 = [0.5; -0.5]; % intial condition for simulation
    intervals = 20; % collect at least as many as the number of unknowns
    iterations = 20; % number of off-policy iterations
    dt = 0.1; % time between samples
elseif system == 3
   [Q, R, K0, A, B, a, f, dim_x, dim_u, K_true, P_true] = params_3d();
    x0 = 10*(rand(3,1)-0.5*ones(3,1)); % intial condition for simulation
    intervals = 30; % collect at least as many as the number of unknowns
    iterations = 20; % number of off-policy iterations
    dt = 0.1; % length of integration interval
else % system == 6
    [Q, R, K0, A, B, a, f, dim_x, dim_u, K_true, P_true] = params_6d(); % get system parameters
    x0 = (rand(6,1)-0.5*ones(6,1)); % intial condition for simulation
    intervals = 70; % collect at least as many as the number of unknowns
    iterations = 30; % number of off-policy iterations
    dt = 0.005; % time between samples

end
%K_true % show true K and P
%P_true

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate the system to collect interval data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nsamples = size(nt,2); % number of different sampling frequencies to try

[kron_x_diffs, x_diffs, xxs, xus, xs, us] = getdata_exact (dim_x, dim_u, intervals, dt, nt(1,end), x0, A, B, K0, a, f);

[rk_adp_exact, K_iterr_exact, P_iterr_exact, K_err_adp_exact, P_err_adp_exact] = adp(dim_x, dim_u, K0, iterations, kron_x_diffs, xxs, xus, Q, R, K_true, P_true);

[rk_id_exact, K_err_id_exact, P_err_id_exact] = sysid(dim_x, dim_u, x_diffs, xs, us, Q, R, K_true, P_true);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate the system to collect sampled data then do quadrature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% arrays for data from ADP
kron_x_diffs_trap = zeros(intervals,dim_x*dim_x,nsamples);
xxs_trap = zeros(intervals,dim_x*dim_x,nsamples);
xus_trap = zeros(intervals,dim_x*dim_u,nsamples);

% arrays for data from system ID
x_diffs_trap = zeros(intervals, dim_x,nsamples);
xs_trap = zeros(intervals,dim_x,nsamples);
us_trap = zeros(intervals,dim_u,nsamples);
   
for j = 1:nsamples % simulate and do quadrature with trapezoid 
    [kron_x_diffs, x_diffs, xxs, xus, xs, us] = getdata_trapezoid (dim_x, dim_u, intervals, dt, nt(1,j), x0, A, B, K0, a, f);
    % save data
    kron_x_diffs_trap(:,:,j) = kron_x_diffs;
    x_diffs_trap(:,:,j) = x_diffs;
    xxs_trap(:,:,j) = xxs;
    xus_trap(:,:,j) = xus;
    xs_trap(:,:,j) = xs;
    us_trap(:,:,j) = us;
end

% arrays for data from ADP
kron_x_diffs_cheb = zeros(intervals,dim_x*dim_x,nsamples);
xxs_cheb = zeros(intervals,dim_x*dim_x,nsamples);
xus_cheb = zeros(intervals,dim_x*dim_u,nsamples);

% arrays for data from system ID
x_diffs_cheb = zeros(intervals, dim_x,nsamples);
xs_cheb = zeros(intervals,dim_x,nsamples);
us_cheb = zeros(intervals,dim_u,nsamples);
   
for j = 1:nsamples  % simulate and do quadrature with chebyshev
    [kron_x_diffs, x_diffs, xxs, xus, xs, us] = getdata_chebyshev (dim_x, dim_u, intervals, dt, nt(1,j), x0, A, B, K0, a, f);
    % save data
    kron_x_diffs_cheb(:,:,j) = kron_x_diffs;
    x_diffs_cheb(:,:,j) = x_diffs;
    xxs_cheb(:,:,j) = xxs;
    xus_cheb(:,:,j) = xus;
    xs_cheb(:,:,j) = xs;
    us_cheb(:,:,j) = us;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do adp on a variable number of samples per interval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[rks_adp_trap, K_errs_adp_trap, P_errs_adp_trap, max_iter_err_trap] = adp_samples(nsamples, iterations, dim_x, dim_u, K0, kron_x_diffs_trap, xxs_trap, xus_trap, Q, R, K_true, P_true); % using trapezoid data
max_iter_err_trap  % show error

[rks_adp_cheb, K_errs_adp_cheb, P_errs_adp_cheb, max_iter_err_cheb] = adp_samples(nsamples, iterations, dim_x, dim_u, K0, kron_x_diffs_cheb, xxs_cheb, xus_cheb, Q, R, K_true, P_true); % using chebyshev data
max_iter_err_cheb  % show error

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do system ID on a variable number of samples per interval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[rks_id_trap, K_errs_id_trap, P_errs_id_trap] = sysid_samples(nsamples, dim_x, dim_u, x_diffs_trap, xs_trap, us_trap, Q, R, K_true, P_true); % using trapezoid data

[rks_id_cheb, K_errs_id_cheb, P_errs_id_cheb] = sysid_samples(nsamples, dim_x, dim_u, x_diffs_cheb, xs_cheb, us_cheb, Q, R, K_true, P_true); % using chebyshev data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tiledlayout(2,1,'TileSpacing','compact','Padding','compact')
nexttile
hold on
h(1) =yline(P_err_adp_exact,'--', 'Color', [0 0.64 0.73],'LineWidth',1);
h(2) =yline(K_err_adp_exact,'-','Color',[0.73 0 0],'LineWidth',1);
for i = 1:nsamples
    if rks_adp_trap(i)
        h(3) = loglog(nt(1,i), P_errs_adp_trap(i), 'diamond', 'MarkerSize',5,'Color',[0 0.49 1],'LineWidth',2);
        h(4) = loglog(nt(1,i), K_errs_adp_trap(i), '+', 'MarkerSize',5,'Color',[0.76 0 0.6],'LineWidth',2);
    end
    if rks_adp_cheb(i)
        h(5) = loglog(nt(1,i), P_errs_adp_cheb(i), 'square', 'MarkerSize',5,'Color',[0.21 0.66 0],'LineWidth',2);
        h(6) = loglog(nt(1,i), K_errs_adp_cheb(i), 'x', 'MarkerSize',5,'Color',[0.77 0.52 0],'LineWidth',2);

    end
end
set(gca,'fontname','times')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlim([nt(1,1), nt(1,end)])
%ylim([1e-4, 1])
%xticks(nt)
set(gca,'Xticklabel',[]) 
set(gca,'fontsize',14)
grid on
ylabel('ADP error')
lgd=legend(h,{'$$P$$', '$$K$$','$$\tilde P^d$$ trap.','$$\tilde K^d$$ trap.','$$\tilde P^d$$ Cheb.','$$\tilde K^d$$ Cheb.'},'location','northeast','interpreter','latex','Orientation','vertical');
lgd.NumColumns = 3;
nexttile
hold on
h(5) =yline(P_err_id_exact,'--', 'Color', [0 0.64 0.73],'LineWidth',1);
h(6) =yline(K_err_id_exact,'-','Color',[0.73 0 0],'LineWidth',1);
for i = 1:nsamples
    if rks_adp_trap(i)
        h(1) = loglog(nt(1,i), P_errs_id_trap(i), 'diamond', 'MarkerSize',5,'Color',[0 0.49 1],'LineWidth',2);
        h(2) = loglog(nt(1,i), K_errs_id_trap(i), '+', 'MarkerSize',5,'Color',[0.76 0 0.6],'LineWidth',2);
    end
    if rks_adp_cheb(i)
        h(3) = loglog(nt(1,i), P_errs_id_cheb(i), 'square', 'MarkerSize',5,'Color',[0.21 0.66 0],'LineWidth',2);
        h(4) = loglog(nt(1,i), K_errs_id_cheb(i), 'x', 'MarkerSize',5,'Color',[0.77 0.52 0],'LineWidth',2);

    end
end

xlim([nt(1,1), nt(1,end)])
%ylim([1e-11, 1e-7])
%xticks([1:intervals])
set(gca,'fontsize',14)
ax = gca;
set(gca,'fontname','times')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
labels = string(ax.XAxis.TickLabels);
ax.XAxis.TickLabels = labels;
grid on
xlabel('Samples per interval')
ylabel('System identification error')
%lgd=legend(h,{'$$\tilde P^d$$ trap.','$$\tilde K^d$$ trap.','$$\tilde P^d$$ Cheb.','$$\tilde K^d$$ Cheb.', '$$P$$', '$$K$$'},'location','northeast','interpreter','latex','Orientation','vertical');
lgd.NumColumns = 3;
set(gcf, 'Position',  [100, 100, 600, 380])