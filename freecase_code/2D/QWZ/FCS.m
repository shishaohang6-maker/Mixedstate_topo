function y=FCS(L1,L2,lambda,beta)
n=2;
N=n*L1;
H=zeros(N,N,L2);%储存哈密顿量
D=zeros(N,L2);%储存本征值
% M=zeros(N,N);
% V=zeros(N);
num=ones(L1,n);
num2=ones(N,2);
nc=0;
for ii=1:L1
    for jj=1:n
        nc=nc+1;
        num(ii,jj)=nc;
        num2(nc,1)=ii;
        num2(nc,2)=jj;
    end
end
for jj=1:L2
    for ii=1:(L1-1)
    n1=num(ii,1);
    n2=num(ii,2);
    nn1=num(ii+1,1);
    nn2=num(ii+1,2);
    H(nn1,n1,jj)=0.5*lambda;
    H(n1,nn1,jj)=0.5*lambda;
    H(nn2,n2,jj)=-0.5*lambda;
    H(n2,nn2,jj)=-0.5*lambda;
    H(nn1,n2,jj)=  sqrt(-1)*0.5*lambda;
    H(n2,nn1,jj)= -sqrt(-1)*0.5*lambda;
    H(nn2,n1,jj)=  sqrt(-1)*0.5*lambda;
    H(n1,nn2,jj)= -sqrt(-1)*0.5*lambda;
    end
end

for jj=1:L2
    k=2*pi*jj/L2;
    for ii=1:L1
        n1=num(ii,1);
        n2=num(ii,2);
        H(n2,n1,jj)= sqrt(-1)*lambda*sin(k);
        H(n1,n2,jj)=-sqrt(-1)*lambda*sin(k);
        H(n1,n1,jj)= 1+lambda*cos(k);
        H(n2,n2,jj)=-1-lambda*cos(k);
    end
end




for jj=1:L2
    M=H(:,:,jj);
    V=eig(M);
    D(:,jj)=sort(V);
end

ky=zeros(L2);
for jj=1:L2
    ky(jj)=2*pi*jj/L1;
end

% figure(1);
% hold on;
% for ii=1:N
%     plot(kx,D(ii,:),'k')
% end
% title(['t2=',num2str(t2),',','m=',num2str(m)])
% hold off;

lin=100;
x=zeros(1,lin);
y=zeros(1,lin);
dy=zeros(1,lin);
ddy=zeros(1,lin);
for nn=1:lin
%    x(nn)=3*pi/4 + (pi/2)*nn/lin;
    x(nn)=2*pi*nn/lin;
    for ii=1:N
        for jj=1:L2
              y(nn)=y(nn)+log(sqrt(1+2*exp(beta*D(ii,jj))*cos(x(nn))+exp(2*beta*D(ii,jj)))/(1+exp(beta*D(ii,jj))))/(L1*L2);
            dy(nn)=dy(nn)-((1/L1)*(1/L2))*exp(beta*D(ii,jj))*sin(x(nn))/(1+2*exp(beta*D(ii,jj))*cos(x(nn))+exp(2*beta*D(ii,jj)));
           ddy(nn)=ddy(nn)-(1/L1)*(1/L2)*exp(beta*D(ii,jj))*((exp(2*beta*D(ii,jj))+1)*cos(x(nn))+2*exp(beta*D(ii,jj)))/(1+2*exp(beta*D(ii,jj))*cos(x(nn))+exp(2*beta*D(ii,jj)))^2;
        end
    end
end

y=ddy(lin/2);

end