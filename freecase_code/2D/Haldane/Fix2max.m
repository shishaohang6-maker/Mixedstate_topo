% === 固定参数 ===
t1 = 1.0;
t2 = 0.1;
m = 0.0;
L2 = 12;
beta_list = [10.0];
L_list = 7:8:87;  % L1 的范围

% === 颜色设置（可按需要扩展）===
colors = ['r', 'b', 'k'];  

% === 准备绘图 ===
figure;
hold on;

for idx = 1:length(beta_list)
    beta = beta_list(idx);
    a = zeros(1, length(L_list));  % F''(\pi)
    b = zeros(1, length(L_list));  % 其他信息（可忽略）

    for n = 1:length(L_list)
        [a(n), b(n)] = Findmax(L_list(n), L2, t1, t2, m, beta);
    end
% === 设定异常值阈值（可调） ===
threshold = 10000.0;   % 比如超过 1000 的值认为是异常

% === 剔除异常值 ===
valid_idx = abs(a) < threshold;  % 逻辑索引：只保留绝对值小于阈值的点
L_valid = L_list(valid_idx);
a_valid = a(valid_idx);

scatter(L_valid, a_valid, 36, colors(idx), 'filled', ...
    'DisplayName', ['$\beta$ = ', num2str(beta)]);
end

% === 图像标注与格式 ===
xlabel('$L_1$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$F''''(\pi)$', 'Interpreter', 'latex', 'FontSize', 14);
title(['$t_1$ = ', num2str(t1), ', $t_2$ = ', num2str(t2), ...
       ', $m$ = ', num2str(m), ', $L_2$ = ', num2str(L2)], ...
       'Interpreter', 'latex');
legend('show', 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'best');
grid on;
hold off;
ammm=zeros(length(L_list),2);
ammm(:,1)=L_list;
ammm(:,2)=a;
xx=ammm(1,1);
yy=ammm(1,2);
ammm(:,1)=ammm(:,1)-xx;
ammm(:,2)=(ammm(:,2)-yy)*beta*beta;
