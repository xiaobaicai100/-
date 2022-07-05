clear all;
clc;
a=1;b=1;
x=-2:0.01:2;
t=[0.02,0.1,0.3];
n=length(x);
wt=zeros(3,n);
for i=1:1:3
wt(i,:)=1/2*(erf((1-x+a*t(i))/sqrt(4*b*t(i)))+erf((1+x-a*t(i))/sqrt(4*b*t(i))));
end
hold on
plot(x,wt(1,:),'g');
plot(x,wt(2,:),'b');
plot(x,wt(3,:),'r');
legend('t=0.02','t=0.1','t=0.3');