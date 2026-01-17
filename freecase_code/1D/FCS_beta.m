L = 32;
n = 2;
t = 1;
dt = 0.1;              % 固定 δt
beta_list = [1.0, 2.0];  % 不同 beta
N = n * L;
lin = 2000;
theta = linspace(0, 2*pi, lin);
colors = lines(length(beta_list));  % 科研风格配色

figure;
hold on;

for b = 1:length(beta_list)
    beta = beta_list(b);
    H = zeros(N, N);
    num = zeros(L, n);
    nc = 0;

    for ii = 1:L
        for jj = 1:n
            nc = nc + 1;
            num(ii, jj) = nc;
        end
    end

    for nn = 1:(L-1)
        n1 = num(nn, 1);
        n2 = num(nn, 2);
        nn1 = num(nn+1, 1);
        H(n1, n2) = t + dt;
        H(n2, n1) = t + dt;
        H(nn1, n2) = t - dt;
        H(n2, nn1) = t - dt;
    end

    H(num(L,1), num(L,2)) = t + dt;
    H(num(L,2), num(L,1)) = t + dt;

    D_o = eig(H);
    D_o = beta * D_o;   % 注意现在是 β * D_o
    y_o = zeros(1, lin);

    for nn = 1:lin
        for ii = 1:N
            y_o(nn) = y_o(nn) + ...
                log( sqrt(1 + 2*exp(-D_o(ii))*cos(theta(nn)) + exp(-2*D_o(ii))) ...
                / (1 + exp(-D_o(ii))) ) / L;
        end
    end

    plot(theta,  y_o, 'LineWidth', 2, ...
         'DisplayName', ['\beta = ', num2str(beta)], ...
         'Color', colors(b,:));
end

hold off;
xlabel('\theta', 'FontSize', 14, 'Interpreter', 'tex');
ylabel('F(\theta)', 'FontSize', 14, 'Interpreter', 'tex');
legend('show', 'FontSize', 12, 'Location', 'best');
set(gca, 'FontSize', 12);
box on;
grid on;
set(gca, 'LineWidth', 1.2);
set(gcf, 'Position', [100, 100, 800, 600]);

