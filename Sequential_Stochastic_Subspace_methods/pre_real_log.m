%读取真实数据集（用于logistics回归） 并进行标准化等处理 此时未分训练集与测试集 n1=a
clc
clear
XXr=xlsread('X2_log_2.xlsx');
XX=XXr';
[a,b]=size(XX);
Xt=[];
Yt=[];
for i=1:a
    Xt(i,:)=XX(i,1:b-1);
    Yt(i,1)=XX(i,b);
end
n1=a;
%%%%%%%将X标准化%%%%%%%%%%
for i=1:b-1
    Xt(:,i)=(Xt(:,i)-mean(Xt(:,i)))/(std(Xt(:,i))*sqrt((a-1)/a));
end
Xt(:,b)=ones(a,1);