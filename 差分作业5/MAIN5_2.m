clear all;
clc;
delt=1/10;%ʱ����
tal=1/10;%x��y�Ŀռ���
t=0:delt:1;
x=-1:tal:1;
y=-1:tal:1;
mu=tal/delt/delt;

U=zeros(length(x),length(y));%��ʼ����
for i=1:length(x)
    for j=1:length(y)
         U(i,j)=sin(x(i)-0.5*y(j))*cosh(x(i)+y(j));
    end
end

for it=2:length(t)  %%%%��ʼѭ��һ��һ��ļ���
    
    U2=zeros(length(x),length(y));%%%�洢n+1/2�����ֵ
    U3=zeros(length(x),length(y));%%%�洢n+1�����
    for j=1:length(y)%%%�߽�����
        U3(1,j)=exp(1.5*(t(1)+(it-1)*delt))*sin(x(1)-0.5*y(j))*cosh(x(1)+y(j));
        U3(length(x),j)=exp(1.5*(t(1)+(it-1)*delt))*sin(x(length(x))-0.5*y(j))*cosh(x(length(x))+y(j));
    end
    for i=2:length(x)-1
        U3(i,1)=exp(1.5*(t(1)+(it-1)*delt))*sin(x(i)-0.5*y(1))*cosh(x(i)+y(1));
        U3(i,length(y))=exp(1.5*(t(1)+(it-1)*delt))*sin(x(i)-0.5*y(length(y)))*cosh(x(i)+y(length(y)));
    end
    
     for j=2:length(y)-1%%%��n+1/2��ı߽�����
         for i=1:length(x)-1:length(x)
             a=(mu+1/6)*(  U(i,j)+1/2*(2*mu+1/6)*(U(i,j+1)-2*U(i,j)+U(i,j-1))     );
             b=(mu-1/6)*(  U3(i,j)-1/2*(2*mu-1/6)*(U3(i,j+1)-2*U3(i,j)+U3(i,j-1))     );
            U2(i,j)=1/(2*mu)*(a+b);
         end
     end
      
        A=ones(1,3);
        A(1)= -1/2*(mu-1/6);      %�洢a
        A(2)= 1+(mu-1/6);  %�洢b
        A(3)=-1/2*(mu-1/6);      %�洢c
        d=zeros(1,length(x)-2);
        for j=2:1:(length(y)-1)
             v=[U2(1,j),U2(length(x),j)]';
             for i=1:length(d)
                 d(i)=U(i+1,j)+1/2*(2*mu+1/6)*( U(i+1,j-1)-2*U(i+1,j)+U(i+1,j+1) );
             end
              [U2(:,j)]=TA(A,v,d); %%v��d�ڲ�ͬ�����в�ͬ����Ҫ�ֱ���
        end
        
         A(1)= -1/2*(2*mu-1/6);      %�洢a
        A(2)= 1+(2*mu-1/6);  %�洢b
        A(3)=-1/2*(2*mu-1/6);      %�洢c
        d1=zeros(1,length(y)-2);
        for i=2:1:(length(x)-1)
           
             v=[U3(i,1),U3(i,length(y))]';
             for j=1:length(d)
                 d1(j)=U2(i,j+1)+1/2*(mu+1/6)*( U2(i+1,j+1)-2*U2(i,j+1)+U2(i-1,j+1) );
             end
              [U3(i,:)]=TA(A,v,d1); %%v��d�ڲ�ͬ�����в�ͬ����Ҫ�ֱ���
        end
    U=U3;
end
Ud=zeros(length(x),length(y));%��ʼ����
for i=1:length(x)
    for j=1:length(y)
         Ud(i,j)=exp(1.5)*sin(x(i)-0.5*y(j))*cosh(x(i)+y(j));
    end
end
Ue=U-Ud;
surf(x,y,Ue);
xlabel('x');
ylabel('y');
zlabel('U');
