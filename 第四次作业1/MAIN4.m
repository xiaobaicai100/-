clear all��
clc
te=1;%�����ʱ��

b=0.1;
tau=10;
y=0;

h=0.1;
nu=1/2;
k=nu*h*h;%%%%�����޸�
x=-1:h:1;
t=0:k:te;
ujingque=exp(-(x-y).^2/(4*b*(te+tau)));%��ȷ��



%%%%%%%%%%%C��N��ʽ
V=zeros(length(t),length(x));
V(1,:)=exp(-(x-y).^2/(4*b*(0+tau)));%��ֵ����
A=ones(1,3);
A(1)= -b*nu/2;       %�洢a
A(2)= 1+b*nu  ;  %�洢b
A(3)=-b*nu/2;      %�洢c
d=zeros(1,length(x)-2);
for i=2:1:length(t)%%��ʼ��ÿ��ʱ������Ӧ��ֵ
    v=[exp(-(-1-y)^2/(4*b*(t(i)+tau))),exp(-(1-y)^2/(4*b*(t(i)+tau)))]';
    for j=1:length(d)
        d(j)=1/2*nu*b*V(i-1,j+2)+(1-b*nu)*V(i-1,j+1)+1/2*nu*b*V(i-1,j);
    end
   [V(i,:)]=TA(A,v,d); %%v��d�ڲ�ͬ�����в�ͬ����Ҫ�ֱ���
end
%%%%%%%%�ŵ��Ը�ʽ
V2=zeros(length(t),length(x));
V2(1,:)=exp(-(x-y).^2/(4*b*(0+tau)));%��ֵ����
for i=2:1:length(t)%%��ʼ��ÿ��ʱ������Ӧ��ֵ
    V2(i,1)=exp(-(-1-y).^2/(4*b*(t(i)+tau)));
    V2(i,length(x))=(exp(-(1-y).^2/(4*b*(t(i)+tau))));
    for j=2:length(x)-1
        V2(i,j)=(1-2*b*nu)*V2(i-1,j)+b*nu*(V2(i-1,j-1)+V2(i-1,j+1));
    end
end
%%%%%%%%%%%back_time central-space��ʽ
V3=zeros(length(t),length(x));
V3(1,:)=exp(-(x-y).^2/(4*b*(0+tau)));%��ֵ����
A=ones(1,3);
A(1)= -b*nu;       %�洢a
A(2)= 1+2*b*nu  ;  %�洢b
A(3)=-b*nu;      %�洢c
d=zeros(1,length(x)-2);
for i=2:1:length(t)%%��ʼ��ÿ��ʱ������Ӧ��ֵ
    v=[exp(-(-1-y)^2/(4*b*(t(i)+tau))),exp(-(1-y)^2/(4*b*(t(i)+tau)))]';
    for j=1:length(d)
        d(j)=V(i-1,j+1);
    end
   [V3(i,:)]=TA(A,v,d); %%v��d�ڲ�ͬ�����в�ͬ����Ҫ�ֱ���
end
%%%%%DFF��ʽ
V4=zeros(length(t),length(x));
V4(1,:)=exp(-(x-y).^2/(4*b*(0+tau)));%��ֵ����
V4(2,:)=exp(-(x-y).^2/(4*b*(k+tau)));%��ֵ����
for i=3:1:length(t)%%��ʼ��ÿ��ʱ������Ӧ��ֵ
    V4(i,1)=exp(-(-1-y)^2/(4*b*(t(i)+tau)));
    V4(i,length(x))=(exp(-(1-y)^2/(4*b*(t(i)+tau))));
    for j=2:length(x)-1
        V4(i,j)=1/(1+2*b*nu)*(2*b*nu*(V4(i-1,j+1)+V4(i-1,j-1))+(1-2*b*nu)*V4(i-2,j));
    end
end

hold on;
plot(x,(V4(length(t),:))','b');
plot(x,(V(length(t),:))','*');
plot(x,ujingque,'r');
