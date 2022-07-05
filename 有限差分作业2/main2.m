clear all
h=1/10;
lamda=1;
k=h*lamda;
x=-1:h:1;
t=0:k:1.2;
V=zeros(length(t),length(x));%第一个变量对应时间的个数，第二个对应空间13*21的矩阵
M=zeros(length(t),length(x)+1);%存储第二种情况下的数值
M(1,:)=[sin(pi*x),sin(pi*(-1+h))];
V(1,:)=sin(pi*x);%初值条件
A=ones(3,length(x)-1);
A(1,:)=-1/4*A(1,:);
A(2,:)=(1-k/2)*A(2,:);
A(3,:)=1/4*A(3,:);%求出将abc放在A中这里实际就是常数
d=zeros(1,length(x)-2);
d1=zeros(1,length(x)-1);
for i=2:1:length(t)%%开始求每个时间所对应的值
    v=[(1+t(i))*sin(pi*(-1-t(i))),(1+t(i))*sin(pi*(1-t(i)))]';
    for j=1:length(d)
       d(j)=-1/4*V(i-1,j+2)+(1+k/2)*V(i-1,j+1)+1/4*V(i-1,j)+k/2*(-t(i)*sin(pi*(x(j+1)-t(i)))-t(i-1)*sin(pi*(x(j+1)-t(i-1)))); 
    end
   [V(i,:)]=TA(A,v,d); %%v，d在不同步骤中不同，需要分别求
  % M(i,:)=PTS(A,d);
end
for i=2:1:length(t)%%开始求每个时间所对应的值
    v=[(1+t(i))*sin(pi*(-1-t(i))),(1+t(i))*sin(pi*(1-t(i)))]';
    for j=1:length(d1)
       d1(j)=-1/4*M(i-1,j+2)+(1+k/2)*M(i-1,j+1)+1/4*M(i-1,j)+k/2*(-t(i)*sin(pi*(x(j+1)-t(i)))-t(i-1)*sin(pi*(x(j+1)-t(i-1)))); 
    end
    M(i,:)=PTS(A,d1);
end
u=2.2*sin(pi*(x-1.2));%精确解

% hold on;
% plot(x,u,'r');
% plot(x,(V(length(t),:))','b');
% plot(x,(M(length(t),1:length(x)))','k');

hold on;
plot(x,(V(length(t),:)-u)','b');
plot(x,(M(length(t),1:length(x))-u)','k');

