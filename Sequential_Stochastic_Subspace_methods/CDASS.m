%CDASS
%%%%%%%%%%%%这是坐标下降的第一种更新方式 naive朴素更新 与m2对比 m2是其全坐标下降方法%%%%%%%%%%%%%%%%
clc %%%%%%%%%%%%%随机子空间的坐标下降 lasso  那些未被选入子空间的坐标不会设为0
clearvars -except Xt Yt n1 a b J0 J1 J2 J3 ws FF FD1 FD3
l=0.1;
e0=10;
e=0.0000005; %最终伊普西龙的值
ksp=50;%每d/ksp次进行一次计数
%%%%%%初始时刻 概率分子分母的值%%%%%%%%%%%%%%
u0=11;
v0=21;
u=u0*ones(1,b-1);
v=v0*ones(1,b-1);

S=zeros(1,b-1);%%%%%子空间坐标索引 1表示选取 0表示不选取
%%%%%%%%%%%上下界阈值%%%%%%%%%%%%
eta1=0.6; %%%真实数据选0.6和0.5
eta2=0.35;
E1=ones(1,b-1)*eta1;
E2=ones(1,b-1)*eta2;

w=zeros(1,b-1); %初值
F=[];
t=1;
F(t)=(Yt-Xt*w')'*(Yt-Xt*w')/(2*n1)+l*sum(abs(w));

%%%%%%%更新方式（1）的所需变量%%%%%%%%
V=Yt-Xt*w';
XX=Xt'*Xt;%12月
VX=Xt'*V;%12月
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kd=0;
kf=1;
sd2=[];
sd2(1)=b-1;
%while (e0>e)
while (F(t)>=FF)
    
    %%%%%%%%%%%12月概率改进式子%%%%%%%%%%%%
    gl=binornd(1,eta2,1,b-1);
    %%%%%%%%%%%12月概率改进式子%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%子空间选择%%%%%%%%%%%%%%%%%
    %S=zeros(1,b-1);
    %for i=1:b-1
    %    if u(i)/v(i)>eta1
    %        S(i)=1;
    %    elseif eta1>=u(i)/v(i) && u(i)/v(i)>eta2
    %        %S(i)=binornd(1,u(i)/v(i));
    %        S(i)=gl(i);
    %    end
    %end
   
    
    %%%%%%%%%%%%%%%%%12月改进的新方法 大大提高计算效率%%%%%%%%%%%%%%%%%%
    q=u./v;
    S=(q>E1);
    S=S+gl.*((q>E2)-ones(1,b-1)+(E1>=q));
    %%%%%%%%%%%%%%%该方法用向量四则运算代替循环%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %%%%%%%%%%%%更新方式（1） begin%%%%%%%%%%%%%%
    for i=1:b-1
        if S(i)==0
            continue
        end
        r=VX(i)/n1+w(i);   %%%%%%%%%%%%%12月将V换为VX
        uu=sign(r)*max(0,abs(r)-l);
        if w(i)~=uu
            VX=VX+XX(:,i)*(w(i)-uu); %%%%%%%%%%%%12月改
            w(i)=uu;
        end
        v(i)=v(i)+1;
        if w(i)~=0
            u(i)=u(i)+1;
        end
        kd=kd+1;
        if mod(kd,ceil((b-1)/ksp))==0;
            FD2(kf)=(Yt-Xt*w')'*(Yt-Xt*w')/(2*n1)+l*sum(abs(w));
            sd2(kf)=sum(w~=0);
            kf=kf+1;
            
        end
    end
    %%%%%%%%%%%%更新方式（1） end%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%概率记录%%%%%%%%%%%%%%%%%%%%%
    %P(t,:)=u./v;
    SS(t,:)=S;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    t=t+1;
    F(t)=(Yt-Xt*w')'*(Yt-Xt*w')/(2*n1)+l*sum(abs(w));
    e0=abs(F(t)-F(t-1));
end
%FF=F(t);

D=[];
D=(Yt-Xt*w')'*(Yt-Xt*w');
p1=sum(D)/n1; %训练集平均相对误差


jg1=sqrt(sum(D))/n1;%Res 二范数平均残差和 ||e||_2/n

o1=0;
o2=0;
%for i=1:b-1
%    o1=o1+abs(ws(i)-w(i));
%    o2=o2+abs(ws(i));
%end
%jg2=o1/o2; %Relerr 参数的二范数和

%jg3=0;
%for i=1:b-1
%    if w(i)==0 & ws(i)==0
%        jg3=jg3+1;
%    elseif w(i)~=0 & ws(i)~=0
%        jg3=jg3+1;
%    end
%end
%jg3=jg3/(b-1); %参数恢复的准确率 即参数选择的准确率

%J3=[];
%J3(1)=t;
%J3(2)=jg1;
%J3(3)=jg2;
%J3(4)=jg3;
%J3=J3';

SSS=sum(SS');
%plot(P(:,1:20))
%plot(P)