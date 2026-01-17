L = 16;
n = 2;
t = 1;
beta = 1.0;                         % 固定 beta
dt_list = [-0.1];  % 不同 dt
N = n * L;
lin = 101;
theta = linspace(0, 2*pi, lin);
colors = lines(length(dt_list));  % 科研风格配色

% figure;
% hold on;

for d = 1:length(dt_list)
    dt = dt_list(d);
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
    D_o = beta * D_o;   % β 固定
    y_o = zeros(1, lin);

    for nn = 1:lin
        for ii = 1:N
            y_o(nn) = y_o(nn) + ...
                log( sqrt(1 + 2*exp(-D_o(ii))*cos(theta(nn)) + exp(-2*D_o(ii))) ...
                / (1 + exp(-D_o(ii))) ) / L;
        end
    end

     plot(theta, y_o, 'LineWidth', 2, ...
          'DisplayName', ['\deltat = ', num2str(dt)], ...
          'Color', colors(d,:));
end

% hold off;
% xlabel('\theta', 'FontSize', 14, 'Interpreter', 'tex');
% ylabel('F(\theta)', 'FontSize', 14, 'Interpreter', 'tex');
% legend('show', 'FontSize', 12, 'Location', 'best');
% set(gca, 'FontSize', 12);
% box on;
% grid on;
% set(gca, 'LineWidth', 1.2);
% set(gcf, 'Position', [100, 100, 800, 600]);
a=zeros(length(theta),2);
a(:,1)=theta';
a(:,2)=y_o';