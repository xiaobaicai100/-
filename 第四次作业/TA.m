%%追赶法 Thomas Algorithm 解三对角矩阵
function [ w] = TA( A,v,d )
%A为一个储存系数矩阵系数的矩阵，A的每一行分别代表a，b,c.
%v是一个二维向量，储存着w0和wm的值
%d代表线性方程组增广矩阵的最后一列
n=length(d);%n存储将要求解未知量的个数m-1个
p=zeros(1,n+1);
q=zeros(1,n+1);
w=zeros(1,n+2);%w0放在第一个位置，wi放在第i+1个位置
w(1)=v(1);
w(n+2)=v(2);
q(1)=v(1);%为p1，q1赋初值
for i=1:1:n %先求出pq的值
    p(i+1)=-1/(A(1)*p(i)+A(2))*A(3);
    q(i+1)=1/(A(1)*p(i)+A(2))*(d(i)-A(1)*q(i));
end
for j=n:-1:1
    w(j+1)=p(j+1)*w(j+2)+q(j+1);
end
end
%%%n个未知量最后输出n+2个值
