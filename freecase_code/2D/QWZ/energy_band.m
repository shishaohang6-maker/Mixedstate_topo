% n=64;
% beta=1.0;
% lambda=0.75;
% L(n)=0;
% L2(n)=0;
% L1=20;
% p(n)=0;
% for ii=1:n
%     L(ii)=9+2*(ii-1);
%     L2(ii)=L(ii);
%     p(ii)=FCS(L1,L2(ii),lambda,beta);
% end
% scatter(L,p);
% title(['\lambda=',num2str(lambda)])

L1=256;L2=32;
lambda=0.75;
n=2;
N=n*L2;
H=zeros(N,N,L1);%储存哈密顿量
D=zeros(N,L1);%储存本征值
% M=zeros(N,N);
% V=zeros(N);
num=ones(L2,n);
num2=ones(N,2);
nc=0;
for jj=1:L2
    for ii=1:n
        nc=nc+1;
        num(jj,ii)=nc;
        num2(nc,1)=jj;
        num2(nc,2)=ii;
    end
end

for ii=1:L1
    for jj=1:(L2-1)

        n1=num(jj,1);
        n2=num(jj,2);
        nn1=num(jj+1,1);
        nn2=num(jj+1,2);

    H(nn1,n1,ii)=0.5*lambda;
    H(n1,nn1,ii)=0.5*lambda;
    H(nn2,n2,ii)=-0.5*lambda;
    H(n2,nn2,ii)=-0.5*lambda;
    H(nn1,n2,ii)=  0.5*lambda;
    H(n2,nn1,ii)=  0.5*lambda;
    H(nn2,n1,ii)= -0.5*lambda;
    H(n1,nn2,ii)= -0.5*lambda;

    end
end

for ii=1:L1
    k=2*pi*ii/L1;
    for jj=1:L2

        n1=num(jj,1);
        n2=num(jj,2);

        H(n2,n1,ii)=lambda*sin(k);
        H(n1,n2,ii)=lambda*sin(k);

        H(n1,n1,ii)= 1+lambda*cos(k);
        H(n2,n2,ii)=-1-lambda*cos(k);

    end
end




for ii=1:L1
    M=H(:,:,ii);
    V=eig(M);
    D(:,ii)=sort(V);
end

kx=zeros(L1);
for ii=1:L1
    kx(ii)=2*pi*ii/L1;
end

figure(1);
hold on;
for ii=1:N
    plot(kx,D(ii,:),'k')
end
hold off;