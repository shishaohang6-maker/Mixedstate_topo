function [p,q]=Findmax(L1,L2,t1,t2,m,beta)
n=2;
N=n*L2;
H=zeros(N,N,L1);%储存哈密顿量
D=zeros(N,L1);%储存本征值
D1=zeros(N,L1);%储存本征值（周期性边界条件）
% M=zeros(N,N);
% V=zeros(N);
num=ones(L2,n);
num2=ones(N,2);
nc=0;
for ii=1:L2
    for jj=1:n
        nc=nc+1;
        num(ii,jj)=nc;
        num2(nc,1)=ii;
        num2(nc,2)=jj;
    end
end
%NN
for nn=1:L1
    for ii=1:(N-1)
        H(ii,ii+1,nn)=H(ii,ii+1,nn)+t1;
        H(ii+1,ii,nn)=H(ii+1,ii,nn)+t1;
    end
end
%NNN
for nn=1:L1
    for ii=1:(L2-1)
        n1=num(ii,1);   n2=num(ii,2);
        nn1=num(ii+1,1);nn2=num(ii+1,2);
        H(n1,nn1,nn)=H(n1,nn1,nn)+1i*t2;H(nn1,n1,nn)=H(nn1,n1,nn)-1i*t2;
        H(n2,nn2,nn)=H(n2,nn2,nn)-1i*t2;H(nn2,n2,nn)=H(nn2,n2,nn)+1i*t2;
    end
end
%chemical
for nn=1:L1
    for ii=1:L2
        n1=num(ii,1);
        n2=num(ii,2);
        H(n1,n1,nn)=H(n1,n1,nn)+m;
        H(n2,n2,nn)=H(n2,n2,nn)-m;
    end
end



%NN 
for nn=1:L1
    k=2*pi*nn/L1;
    for ii=1:L2
        n1=num(ii,1);n2=num(ii,2);  
        H(n1,n2,nn)=H(n1,n2,nn)+t1*exp( 1i*k);
        H(n2,n1,nn)=H(n2,n1,nn)+t1*exp(-1i*k);
    end
end
%NNN
for nn=1:L1
    k=2*pi*nn/L1;
    for ii=1:L2
        n1=num(ii,1);n2=num(ii,2);
        H(n1,n1,nn)=H(n1,n1,nn)-t2*2*sin(k);
        H(n2,n2,nn)=H(n2,n2,nn)+t2*2*sin(k);
    end
end


for nn=1:L1
    k=2*pi*nn/L1;
    for ii=1:(L2-1)
        n1=num(ii,1);   n2=num(ii,2);
        nn1=num(ii+1,1);nn2=num(ii+1,2);
        H(n1,nn1,nn)=H(n1,nn1,nn)-1i*t2*exp( 1i*k);
        H(nn1,n1,nn)=H(nn1,n1,nn)+1i*t2*exp(-1i*k);
        H(n2,nn2,nn)=H(n2,nn2,nn)+1i*t2*exp( 1i*k);
        H(nn2,n2,nn)=H(nn2,n2,nn)-1i*t2*exp(-1i*k);
    end
end

for ii=1:L1
    M=H(:,:,ii);
    V=eig(M);
    D(:,ii)=sort(V);
end

% kx=zeros(L1);
% for nn=1:L1
%     kx(nn)=2*pi*nn/L1;
% end
for ii=1:L2
    for jj=1:L1
        k1=2*pi*jj/L1;
        k2=2*pi*ii/L2;
        D1(ii,jj)=sqrt((m-2*t2*(sin(k1)-sin(k2)+sin(k2-k1)))^2+(t1^2)*(1+cos(k1)+cos(k2))^2+(t1^2)*(sin(k1)+sin(k2))^2);
        D1(N-ii+1,jj)=-D1(ii,jj);
    end
end
%figure(1);
% hold on;
% for ii=1:N
%     plot(kx,D(ii,:),'k')
% end
% title(['t2=',num2str(t2),',','m=',num2str(m)])
% hold off;

lin=2;
x=zeros(1,lin);
y=zeros(1,lin);
y1=zeros(1,lin);
 dy=zeros(1,lin);
dy1=zeros(1,lin);
 ddy=zeros(1,lin);
ddy1=zeros(1,lin);
for nn=1:lin
%   x(nn)=3*pi/4 + (pi/2)*nn/lin;
    x(nn)=2*pi*nn/lin;
    for ii=1:N
        for jj=1:L1
              y(nn)=y(nn)+log(sqrt(1+2*exp(beta*D(ii,jj))*cos(x(nn))+exp(2*beta*D(ii,jj)))/(1+exp(beta*D(ii,jj))))/(L1*L2);
             y1(nn)=y1(nn)+log(sqrt(1+2*exp(beta*D1(ii,jj))*cos(x(nn))+exp(2*beta*D1(ii,jj)))/(1+exp(beta*D1(ii,jj))))/(L1*L2);
            dy(nn)=dy(nn)-((1/L1)*(1/L2))*exp(beta*D(ii,jj))*sin(x(nn))/(1+2*exp(beta*D(ii,jj))*cos(x(nn))+exp(2*beta*D(ii,jj)));
           dy1(nn)=dy1(nn)-((1/L1)*(1/L2))*exp(beta*D1(ii,jj))*sin(x(nn))/(1+2*exp(beta*D1(ii,jj))*cos(x(nn))+exp(2*beta*D1(ii,jj)));
           ddy(nn)=ddy(nn)-(1/L1)*(1/L2)*exp(beta*D(ii,jj))*((exp(2*beta*D(ii,jj))+1)*cos(x(nn))+2*exp(beta*D(ii,jj)))/(1+2*exp(beta*D(ii,jj))*cos(x(nn))+exp(2*beta*D(ii,jj)))^2;
          ddy1(nn)=ddy1(nn)-(1/L1)*(1/L2)*exp(beta*D1(ii,jj))*((exp(2*beta*D1(ii,jj))+1)*cos(x(nn))+2*exp(beta*D1(ii,jj)))/(1+2*exp(beta*D1(ii,jj))*cos(x(nn))+exp(2*beta*D1(ii,jj)))^2;
        end
    end
end
%figure(2);
% hold on;
% plot(x,y)
g=ddy;
[p,q]=max(g);
%q=3*pi/4 + (pi/2)*q/lin;
q=2*pi*q/lin;
%hold off;
end