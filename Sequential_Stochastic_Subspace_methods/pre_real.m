%读取真实数据集 并进行标准化等处理 此时未分训练集与测试集 n1=a
clc
clear
XXr=xlsread('X1_lin.xlsx');
XX=XXr';
[a,b]=size(XX);
Xt=[];
Yt=[];
for i=1:a
    Xt(i,:)=XX(i,1:b-1);
    Yt(i,1)=XX(i,b);
end
n1=a;
%%%%%%%将Y中心化 X标准化%%%%%%%%%%
Yt=Yt-mean(Yt);
for i=1:b-1
    Xt(:,i)=(Xt(:,i)-mean(Xt(:,i)))/(std(Xt(:,i))*sqrt((a-1)/a));
end