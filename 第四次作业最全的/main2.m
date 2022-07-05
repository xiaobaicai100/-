clear all；
clc
te=1;%所求的时刻

b=0.1;
tau=10;
y=0;

h=0.1;
nu=1/2;
k=nu*h*h;%%%%可以修改
x=-2:h:2;
t=0:k:te;
ujingque=1/2*(erf((1-x)/sqrt(4*b*te))+erf((1+x)/sqrt(4*b*te)));%精确解
%%%%%初值条件
u0=zeros(1,length(x));
for i=(1/h+1):1:(3/h+1)
    u0(i)=1;
end

%%%%%%%%%%%C—N格式
V=zeros(length(t),length(x));
V(1,:)=u0;%初值条件
A=ones(1,3);
A(1)= -b*nu/2;       %存储a
A(2)= 1+b*nu  ;  %存储b
A(3)=-b*nu/2;      %存储c
d=zeros(1,length(x)-2);
for i=2:1:length(t)%%开始求每个时间所对应的值
    v=[1/2*(erf((1-x)/sqrt(4*b*t(i)))+erf((1+x)/sqrt(4*b*t(i)))),1/2*(erf((1-x)/sqrt(4*b*t(i)))+erf((1+x)/sqrt(4*b*t(i))))]';
    for j=1:length(d)
        d(j)=1/2*nu*b*V(i-1,j+2)+(1-b*nu)*V(i-1,j+1)+1/2*nu*b*V(i-1,j);
    end
   [V(i,:)]=TA(A,v,d); %%v，d在不同步骤中不同，需要分别求
end
%%%%%%%%古典显格式
V2=zeros(length(t),length(x));
V2(1,:)=u0;%初值条件
for i=2:1:length(t)%%开始求每个时间所对应的值
    V2(i,1)=1/2*(erf((1-(-2))/sqrt(4*b*t(i)))+erf((1+(-2))/sqrt(4*b*t(i))));
    V2(i,length(x))=1/2*(erf((1-2)/sqrt(4*b*t(i)))+erf((1+2)/sqrt(4*b*t(i))));
    for j=2:length(x)-1
        V2(i,j)=(1-2*b*nu)*V2(i-1,j)+b*nu*(V2(i-1,j-1)+V2(i-1,j+1));
    end
end
%%%%%%%%%%%back_time central-space格式
V3=zeros(length(t),length(x));
V3(1,:)=u0;%初值条件
A=ones(1,3);
A(1)= -b*nu;       %存储a
A(2)= 1+2*b*nu  ;  %存储b
A(3)=-b*nu;      %存储c
d=zeros(1,length(x)-2);
for i=2:1:length(t)%%开始求每个时间所对应的值
    v=[1/2*(erf((1-x)/sqrt(4*b*t(i)))+erf((1+x)/sqrt(4*b*t(i)))),1/2*(erf((1-x)/sqrt(4*b*t(i)))+erf((1+x)/sqrt(4*b*t(i))))]';
    for j=1:length(d)
        d(j)=V(i-1,j+1);
    end
   [V3(i,:)]=TA(A,v,d); %%v，d在不同步骤中不同，需要分别求
end

%%%%%DFF格式
V4=zeros(length(t),length(x));
V4(1,:)=u0;%初值条件
V4(2,:)=1/2*(erf((1-x)/sqrt(4*b*t(2)))+erf((1+x)/sqrt(4*b*t(2))));%初值条件
for i=3:1:length(t)%%开始求每个时间所对应的值
   V2(i,1)=1/2*(erf((1-(-2))/sqrt(4*b*t(i)))+erf((1+(-2))/sqrt(4*b*t(i))));
    V2(i,length(x))=1/2*(erf((1-2)/sqrt(4*b*t(i)))+erf((1+2)/sqrt(4*b*t(i))));
    for j=2:length(x)-1
        V4(i,j)=1/(1+2*b*nu)*(2*b*nu*(V4(i-1,j+1)+V4(i-1,j-1))+(1-2*b*nu)*V4(i-2,j));
    end
end

hold on;

plot(x,(V(length(t),:))','m');
 plot(x,(V2(length(t),:))','b');
plot(x,(V3(length(t),:))','y');
 plot(x,(V4(length(t),:))','g');
plot(x,ujingque,'r');
legend('C-Nshame','forward-time central-space shame','back_time central-space shame','Du Fort-Frankelshame','exact solution');