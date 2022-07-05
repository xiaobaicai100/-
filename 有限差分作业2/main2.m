clear all
h=1/10;
lamda=1;
k=h*lamda;
x=-1:h:1;
t=0:k:1.2;
V=zeros(length(t),length(x));%��һ��������Ӧʱ��ĸ������ڶ�����Ӧ�ռ�13*21�ľ���
M=zeros(length(t),length(x)+1);%�洢�ڶ�������µ���ֵ
M(1,:)=[sin(pi*x),sin(pi*(-1+h))];
V(1,:)=sin(pi*x);%��ֵ����
A=ones(3,length(x)-1);
A(1,:)=-1/4*A(1,:);
A(2,:)=(1-k/2)*A(2,:);
A(3,:)=1/4*A(3,:);%�����abc����A������ʵ�ʾ��ǳ���
d=zeros(1,length(x)-2);
d1=zeros(1,length(x)-1);
for i=2:1:length(t)%%��ʼ��ÿ��ʱ������Ӧ��ֵ
    v=[(1+t(i))*sin(pi*(-1-t(i))),(1+t(i))*sin(pi*(1-t(i)))]';
    for j=1:length(d)
       d(j)=-1/4*V(i-1,j+2)+(1+k/2)*V(i-1,j+1)+1/4*V(i-1,j)+k/2*(-t(i)*sin(pi*(x(j+1)-t(i)))-t(i-1)*sin(pi*(x(j+1)-t(i-1)))); 
    end
   [V(i,:)]=TA(A,v,d); %%v��d�ڲ�ͬ�����в�ͬ����Ҫ�ֱ���
  % M(i,:)=PTS(A,d);
end
for i=2:1:length(t)%%��ʼ��ÿ��ʱ������Ӧ��ֵ
    v=[(1+t(i))*sin(pi*(-1-t(i))),(1+t(i))*sin(pi*(1-t(i)))]';
    for j=1:length(d1)
       d1(j)=-1/4*M(i-1,j+2)+(1+k/2)*M(i-1,j+1)+1/4*M(i-1,j)+k/2*(-t(i)*sin(pi*(x(j+1)-t(i)))-t(i-1)*sin(pi*(x(j+1)-t(i-1)))); 
    end
    M(i,:)=PTS(A,d1);
end
u=2.2*sin(pi*(x-1.2));%��ȷ��

% hold on;
% plot(x,u,'r');
% plot(x,(V(length(t),:))','b');
% plot(x,(M(length(t),1:length(x)))','k');

hold on;
plot(x,(V(length(t),:)-u)','b');
plot(x,(M(length(t),1:length(x))-u)','k');

