clear; 

% ==== 固定参数 ====
L2 = 16;
lambda = 0.75;
L1_list = 19:8:189;
beta_list = [10.0];


am=zeros(length(L1_list),2);

% figure('Color','w'); hold on;

colors = lines(length(beta_list));   % 科研风配色

for ib = 1:length(beta_list)
    beta = beta_list(ib);

    % 计算 FCSopen_y
    ddy = arrayfun(@(L1) FCSopen_y(L1, L2, lambda, beta), L1_list);

    % 减去起点
    dL1 = L1_list - L1_list(1);
    ddy_shift = (ddy - ddy(1)) * beta^2;

     % 绘制空心散点
     plot(dL1, ddy_shift, 'o', ...
         'MarkerSize', 7, ...
         'MarkerEdgeColor', colors(ib,:), ...
         'MarkerFaceColor', 'none', ...
         'LineWidth', 1.2);
% 
%     % 虚线连接
%     plot(dL1, ddy_shift, '--', ...
%         'Color', colors(ib,:), ...
%         'LineWidth', 1.2);

    % legend 字符串
%     legends{ib} = sprintf('\\beta = %.1f', beta);
end

% % ==== 图形美化 ====
% xlabel('\Delta L_x', 'Interpreter','latex', 'FontSize', 20);
% ylabel('$\mathcal{F}$', 'Interpreter','latex', 'FontSize', 22);
% 
% legend(legends, 'Interpreter','latex', ...
%        'FontSize', 16, 'Location','northwest', ...
%        'Box','on');
% 
% title('(b)', 'FontSize', 26, 'FontWeight','bold');
% 
% set(gca, 'LineWidth', 1.2, 'FontSize', 16);

% 去除负数轴（很重要）
% xlim([0, max(dL1)]);
% ylim([0, max(ddy_shift)*1.05]);
% 
% box on;
am(:,1) = dL1;
am(:,2) = ddy_shift;
