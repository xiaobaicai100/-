%forward-time central-space scheme 
clear all;
h=1/40;%距离间隔
lam=0.8;
k=h*lam;%时间间隔
t=0:k:2.4;
x=-1:h:3;
lt=length(t);
lx=length(x);
v(lx,lt)=0;
for i=1:1:lx  %赋初值
    if((-1+(i-1)*h<(-1/2))||(-1+(i-1)*h>(1/2)))
        v(i,1)=0;
    else
        v(i,1)=cos((-1+(i-1)*h)*pi)*cos((-1+(i-1)*h)*pi);
    end
end
for i=1:lt%为边界赋值
   v(1,i)=0 ;
end
for i=2:lt
 for j=2:lx-1
    v(j,i)=v(j,i-1)+(lam/2)*v(j-1,i-1)-(lam/2)*v(j+1,i-1);
 end
 v(lx,i)=v(lx-1,i);
end
plot(x,(v(:,lt))');