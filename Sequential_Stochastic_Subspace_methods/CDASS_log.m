%%%%%%%%%%%%%%%子空间坐标下降CDASS 新实验 logistics回归
%%%%%%%%%%%%这是logistics回归的随机子空间坐标下降  与m_logistic对比 m_logistic是其全坐标下降方法%%%%%%%%%%%%%%%%
clc %%%%%%%%%%%%%随机子空间的坐标下降 lasso  那些未被选入子空间的坐标不会设为0
clearvars -except Xt Yt n1 a b J0 J1 J2 J3 ws FF FD1 FD3 sd1 sd3
l=0.1;
e0=10;
e=0.000001; %最终伊普西龙的值
ksp=50;%每d/ksp次进行一次计数
%%%%%%初始时刻 概率分子分母的值%%%%%%%%%%%%%%
u0=11;
v0=21;
u=u0*ones(1,b-1);
v=v0*ones(1,b-1);
na=ones(1,a)/a;
%%%%%%%%%%%注意 此时含有常数项 第b个分量为常数项对应的系数
w=zeros(1,b); %初值 
F=[];
t=1;
kkp=Yt*ones(1,b).*Xt;
F(t)=na*log(1+exp(-kkp*w'))+l*sum(abs(w(1:b-1)));


S=zeros(1,b-1);%%%%%子空间坐标索引 1表示选取 0表示不选取
%%%%%%%%%%%上下界阈值%%%%%%%%%%%%
eta1=0.6; %%%真实数据选0.65和0.5
eta2=0.35;
E1=ones(1,b-1)*eta1;
E2=ones(1,b-1)*eta2;

kd=0;
kf=1;
sd2=[];
sd2(1)=b-1;

while (F(t)>=FF)
    %%%%%%%%%%%12月概率改进式子%%%%%%%%%%%%
    gl=binornd(1,eta2,1,b-1);
    %%%%%%%%%%%12月概率改进式子%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%12月改进的新方法 大大提高计算效率%%%%%%%%%%%%%%%%%%
    q=u./v;
    S=(q>E1);
    S=S+gl.*((q>E2)-ones(1,b-1)+(E1>=q));
    %%%%%%%%%%%%%%%该方法用向量四则运算代替循环%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i=1:b-1
        if S(i)==0
            continue
        end
        pj=-kkp(:,i)'*(1./(1+exp(kkp*w')))/a;
        r=w(i)-4*pj;
        w(i)=sign(r)*max(0,abs(r)-4*l);
        
        kd=kd+1;
        if mod(kd,ceil((b-1)/ksp))==0;
            FD2(kf)=na*log(1+exp(-kkp*w'))+l*sum(abs(w(1:b-1)));
            sd2(kf)=sum(w~=0);
            kf=kf+1;
            
        end
        
        v(i)=v(i)+1;
        if w(i)~=0
            u(i)=u(i)+1;
        end
        %%%%%%%log的加速计算%%%%%%%%
        %wt(i)=sign(r)*max(0,abs(r)-l);
        %if wt(i)~=w(i)
        %    kp1=exp(-Xt*w');
        %    pj=kkp+(1./(1+kp1))'*Xt/a;
        %    w(i)=wt(i);
        %end
        %%%%%%%%%%%%%概率更新（新方法中对应的概率）%%%%%%%%%%%%%
    %v(i)=v(i)+1;
     %   if w(i)~=0
     %       u(i)=u(i)+1;
     %   end
    %%%%%%%%%%%%%%%%%%%%%%%%%
    end
    pj=-kkp(:,b)'*(1./(1+exp(kkp*w')))/a;
    w(b)=w(b)-4*pj;
    
    %%%%%%%%%%%%更新方式（2） end%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%更新方式（3）原始梯度下降 begin%%%%%%%%%%%%%%
    %for i=1:b-1
    %    r=Xt(:,i)'*(Yt-Xt*w')/n1+w(i);
    %    w(i)=sign(r)*max(0,abs(r)-l);
    %end
    %%%%%%%%%%%%更新方式（2） end%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%概率记录%%%%%%%%%%%%%%%%%%%%%
    P(t,:)=u./v;
    SS(t,:)=S;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t=t+1
    F(t)=na*log(1+exp(-kkp*w'))+l*sum(abs(w(1:b-1)));
    %e0=abs(F(t)-F(t-1));
end


D=[];

PP1=sign(Yt.*(Xt*w'));
PP2=1./(1+exp(-kkp*w'));
%SSS=sum(SS');
%plot(P(:,1:20))
%plot(P)