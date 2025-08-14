clear all
rng(6,"twister");
warning('off')
format long 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% choose system 2 or 3 or 6 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

system = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters for simulation and data collection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if system == 2
    [Q, R, K0, A, B, a, f, dim_x, dim_u, K_true, P_true] = params_2d();
    x0 = [0.5; -0.5]; % intial condition for simulation
    intervals = 30; % collect at least as many as the number of unknowns
    iterations = 20; % number of off-policy iterations
    dt = 0.1; % time between samples
    nt = 100; % number of samples per interval
elseif system == 3
   [Q, R, K0, A, B, a, f, dim_x, dim_u, K_true, P_true] = params_3d();
    x0 = 10*(rand(3,1)-0.5*ones(3,1)); % intial condition for simulation
    intervals = 30; % collect at least as many as the number of unknowns
    iterations = 20; % number of off-policy iterations
    dt = 0.1; % length of integration interval
    nt = 100; % number of samples per interval
else % system == 6
    [Q, R, K0, A, B, a, f, dim_x, dim_u, K_true, P_true] = params_6d(); % get system parameters
    x0 = (rand(6,1)-0.5*ones(6,1)); % intial condition for simulation
    intervals = 70; % collect at least as many as the number of unknowns
    iterations = 30; % number of off-policy iterations
    dt = 0.005; % time between samples
    nt = 100; % number of samples per interval

end
%K_true % show true K and P
%P_true

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate the system to collect data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[kron_x_diffs, x_diffs, xxs, xus, xs, us] = getdata_exact (dim_x, dim_u, intervals, dt, nt, x0, A, B, K0, a, f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do adp on a variable number of intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[rks_adp, K_errs_adp, P_errs_adp, max_iter_err] = adp_intervals(intervals, iterations, dim_x, dim_u, K0, kron_x_diffs, xxs, xus, Q, R, K_true, P_true);
max_iter_err  % show error

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do system ID on a variable number of intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[rks_id, K_errs_id, P_errs_id] = sysid_intervals(intervals, dim_x, dim_u, x_diffs, xs, us, Q, R, K_true, P_true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tiledlayout(2,1,'TileSpacing','compact','Padding','compact')
nexttile
hold on
for i = 1:intervals
    if rks_adp(i)
        h(1) = semilogy(i, P_errs_adp(i), 'o', 'MarkerSize',5,'Color',[0 0.64 0.73],'LineWidth',2);
        h(2) = semilogy(i, K_errs_adp(i), '.', 'MarkerSize',10,'Color',[0.73 0 0],'LineWidth',2);
    end
end
set(gca,'fontname','times')
set(gca, 'YScale', 'log')
xlim([1 intervals])
legend(h,{'$$P$$','$$K$$'},'location','northeast','interpreter','latex')
set(gca,'Xticklabel',[]) 
set(gca,'fontsize',14)
grid on
ylabel('ADP error')
lgd=legend(h,{'$$P$$', '$$K$$'},'location','northeast','interpreter','latex','Orientation','vertical');
lgd.NumColumns = 1;
nexttile
hold on
for i = 1:intervals

    if rks_id(i)
        h(1) = semilogy(i, P_errs_id(i), 'o', 'MarkerSize',5,'Color',[0 0.64 0.73],'LineWidth',2);
        h(2) = semilogy(i, K_errs_id(i), '.', 'MarkerSize',10,'Color',[0.73 0 0],'LineWidth',2);
       
    end
end
xlim([1 intervals])
set(gca,'fontsize',14)
ax = gca;
set(gca,'fontname','times')
set(gca, 'YScale', 'log')
labels = string(ax.XAxis.TickLabels);
ax.XAxis.TickLabels = labels;
grid on
xlabel('Number of intervals')
ylabel('System identification error')
set(gcf, 'Position',  [100, 100, 600, 380])
