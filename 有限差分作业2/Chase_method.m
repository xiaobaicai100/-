function [ x ] = Chase_method( A, b )
%Chase method ׷�Ϸ������ԽǾ���Ľ�
%   AΪ���ԽǾ����ϵ����bΪ��ʽ�Ҷ˵ĳ��������ֵx��Ϊ���յĽ�
%   ע��A����Ϊ����bһ��ҪΪ������
%% ��׷�Ϸ�����L��U
T = A;
for i = 2 : size(T,1)
    T(i,i-1) = T(i,i-1)/T(i-1,i-1);
    T(i,i) = T(i,i) - T(i-1,i) * T(i,i-1);
end
L = zeros(size(T));
L(logical(eye(size(T)))) = 1;   %�Խ��߸�ֵ1
for i = 2:size(T,1)
    for j = i-1:size(T,1)
        L(i,j) = T(i,j);
        break;
    end
end
U = zeros(size(T));
U(logical(eye(size(T)))) = T(logical(eye(size(T))));
for i = 1:size(T,1)
    for j = i+1:size(T,1)
        U(i,j) = T(i,j);
        break;
    end
end
%% ����matlab����󷽳̵ı���ֱ�����
y = L\b;
x = U\y;
end