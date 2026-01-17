% === 参数设置 ===
L1 = 36;
L2 = 36;
t1 = 1;
t2 = 0.1;
beta = 10.0;
ndis = 1000;

n = 2;
N = n * L2;

% === 格点编号预处理 ===
num = ones(L2, n);
num2 = ones(N, 2);

nc = 0;
for ii = 1:L2
    for jj = 1:n
        nc = nc + 1;
        num(ii, jj) = nc;
        num2(nc, :) = [ii, jj];
    end
end

% === 多个 m 的取值 ===
m_list = [0.6];
colors = ['r', 'b', 'g', 'k', 'm', 'c'];  % 用于区分颜色，按需扩展

% figure;
%hold on;

for m_idx = 1:length(m_list)
    m = m_list(m_idx);
    
    x = zeros(1, ndis);
    y = zeros(1, ndis);

    for nt = 1:ndis

        H = zeros(N, N, L1);
        D = zeros(N, L1);

        % === 构造 NN 项 ===
        for nn = 1:L1
            for ii = 1:(N - 1)
                H(ii, ii+1, nn) = H(ii, ii+1, nn) + t1;
                H(ii+1, ii, nn) = H(ii+1, ii, nn) + t1;
            end
        end

        % === 构造 NNN 项（y方向）===
        for nn = 1:L1
            for ii = 1:(L2 - 1)
                n1 = num(ii,1); n2 = num(ii,2);
                nn1 = num(ii+1,1); nn2 = num(ii+1,2);
                H(n1, nn1, nn) = H(n1, nn1, nn) + 1i * t2;
                H(nn1, n1, nn) = H(nn1, n1, nn) - 1i * t2;
                H(n2, nn2, nn) = H(n2, nn2, nn) - 1i * t2;
                H(nn2, n2, nn) = H(nn2, n2, nn) + 1i * t2;
            end
        end

        % === 化学势项 ===
        for nn = 1:L1
            for ii = 1:L2
                n1 = num(ii,1);
                n2 = num(ii,2);
                H(n1,n1,nn) = H(n1,n1,nn) + m;
                H(n2,n2,nn) = H(n2,n2,nn) - m;
            end
        end

        % === 动量方向上的 NN 项 ===
        for nn = 1:L1
            k = 2 * pi * nn / L1 - 2 * pi * nt / ndis / L1;
            for ii = 1:L2
                n1 = num(ii,1); n2 = num(ii,2);
                H(n1,n2,nn) = H(n1,n2,nn) + t1 * exp(1i * k);
                H(n2,n1,nn) = H(n2,n1,nn) + t1 * exp(-1i * k);
            end
        end

        % === 动量方向上的 NNN 项 ===
        for nn = 1:L1
            k = 2 * pi * nn / L1 - 2 * pi * nt / ndis / L1;
            for ii = 1:L2
                n1 = num(ii,1); n2 = num(ii,2);
                H(n1,n1,nn) = H(n1,n1,nn) - 2 * t2 * sin(k);
                H(n2,n2,nn) = H(n2,n2,nn) + 2 * t2 * sin(k);
            end
        end

        for nn = 1:L1
            k = 2 * pi * nn / L1 - 2 * pi * nt / ndis / L1;
            for ii = 1:(L2 - 1)
                n1 = num(ii,1); n2 = num(ii,2);
                nn1 = num(ii+1,1); nn2 = num(ii+1,2);
                H(n1,nn1,nn) = H(n1,nn1,nn) - 1i * t2 * exp(1i * k);
                H(nn1,n1,nn) = H(nn1,n1,nn) + 1i * t2 * exp(-1i * k);
                H(n2,nn2,nn) = H(n2,nn2,nn) + 1i * t2 * exp(1i * k);
                H(nn2,n2,nn) = H(nn2,n2,nn) - 1i * t2 * exp(-1i * k);
            end
        end

        % === 对每个 k 计算本征值 ===
        for ii = 1:L1
            M = H(:,:,ii);
            V = eig(M);
            D(:,ii) = sort(real(V));
        end

        % === 计算 F 值 ===
        x(nt) = nt / ndis;
        for ii = 1:N
            for jj = 1:L1
                e = exp(beta * D(ii,jj));
                nume = sqrt(1 - 2 * e + e^2);
                deno = 1 + e;
                y(nt) = y(nt) + log(nume / deno) / (L1 * L2);
            end
        end
    end

%     % === 绘图 ===
     plot(x, y, 'Color', colors(m_idx), 'LineWidth', 2, ...
         'DisplayName', sprintf('m = %.2f', m));
end

% xlabel('$\tau$', 'Interpreter', 'latex', 'FontSize', 14);
% ylabel('$F(\pi)$', 'Interpreter', 'latex', 'FontSize', 14);
% legend('show', 'FontSize', 12);
% grid on;

a=zeros(ndis,2);
a(:,1)=x';
a(:,2)=y';
