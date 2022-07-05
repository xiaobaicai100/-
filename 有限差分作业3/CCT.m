a=1;b=1;
x=-2:0.01:2;
t=0.02;
wt1=1/2*(erf((1-x+a*t)/sqrt(4*b*t))+erf((1+x-a*t)/sqrt(4*b*t)));
t=0.1;
wt2=1/2*(erf((1-x+a*t)/sqrt(4*b*t))+erf((1+x-a*t)/sqrt(4*b*t)));
t=0.5;
wt3=1/2*(erf((1-x+a*t)/sqrt(4*b*t))+erf((1+x-a*t)/sqrt(4*b*t)));
hold on
plot(x,wt1,'k');
plot(x,wt2,'b');
plot(x,wt3,'r');

