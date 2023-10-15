clc;
clear;
close all;
global costp Ppv Pwt %costp可能代表成本参数,Ppv和Pwt表示光伏和风力发电功率
%% 算法参数
parameter;                %parameter关键字用于声明函数的参数,使参数名称对函数体可见,这样在函数体内就可以直接使用这些参数变量。
nVar=4*24;                % Number of Decision Variables  用于表示优化问题的决策变量个数
VarMin=[ones(1,24)*Pmt_min, ones(1,24)*Pfc_min, ones(1,24)*Px_min, ones(1,24)*Pb_min];%优化变量的下界限制
VarMax=[ones(1,24)*Pmt_max, ones(1,24)*Pfc_max, ones(1,24)*Px_max, ones(1,24)*Pb_max];%优化变量的上界限制
MaxIt=100;                % Maximum Number of Iterations 最大迭代次数
nPop=500;                 % Population Size (Swarm Size) 进化群体规模

%% 计算
[ bestPosition, fitValue ] =SSA(@objective,nVar,VarMin,VarMax,MaxIt,nPop);
x=bestPosition;
Pmt = x(1:24);            % 燃气轮机功率
Pfc = x(25:48);           % 燃料电池功率
Px = x(49:72);            % 可卸负荷               
Pb = x(73:96);            % 蓄电池功率
t=1:24;

%% 输出光伏出力预测
figure
plot(t,Ppv,'-')
title('光伏发电曲线');
xlabel('时间/小时')
ylabel('功率/kw')
%% 输出风力发电出力预测
figure
plot(t,Pwt,'-')
title('风力发电曲线');
xlabel('时间/小时')
ylabel('功率/kw')
%% 输出电平衡结果
figure
hold on 
Pb_po=max(Pb,0);
Pb_ne=min(Pb,0);
positive=[Pmt', Pfc', Pb_po',Px',Ppv',Pwt'];
negative=[ Pb_ne'];
bar(positive,'stack');
bar(negative,'stack');
plot(t, Pl, 'ok-');
title('电平衡');
legend('Pmt燃机','Pfc燃电池','Pbdis电池放电','Px可去负荷','Ppv光伏发电','Pwt风电' ,'Pbch电池充电','Pl总负荷');
grid on
hold off
xlabel('时间/小时')
ylabel('功率/kw')


figure
plot(t, Pl, 'ok-');
set(gca,'xtick',1:24)
title('总负荷');
grid on
xlabel('时间/小时')
ylabel('功率/kw')
%% 输出各部分出力结果
figure
plot(t,Pmt,'ok-')
hold on
plot(t,Pfc,'-*')
hold on
plot(t,Px,'-')
hold on
plot(t,Pb,'-.')
legend('Pmt燃机出力','Pfc燃电池出力','Px可去负荷','Pb电池出力');
title('出力图');
xlabel('时间/小时')
ylabel('功率/kw')
%% 计算每小时运行费用
% 预分配
eta_mt = zeros(1,24);
eta_fc = zeros(1,24);
Pmth = zeros(1,24);
Umt = zeros(1,24);
Ufc = zeros(1,24);


%% 运行模型
for t=1:24
    % 燃气轮机热功率
    %下面这个是微型燃气轮机效率计算公式
    eta_mt(t) = 0.0753*(Pmt(t)/65)^3 - 0.3095*(Pmt(t)/65)^2 + 0.4174*(Pmt(t)/65) + 0.1068;
    Pmth(t) = ((Pmt(t)*(1-eta_mt(t)-eta_l))/eta_mt(t))*eta_h*Coph;
     % 燃料电池功率
    eta_fc(t) = -0.0023*Pfc(t) + 0.674;
end



%% 启停 
for t=1:24
    if Pmt(t)>0
        Umt(t) = 1;
    end
   
    if Pfc(t)>0
        Ufc(t) = 1;
    end
end

for t=1:24 
    if t==1
        Cst(t) = Cst_mt*max(0, Umt(t)-Uinit)  + Cst_fc*max(0, Ufc(t)-Uinit);
    else
        Cst(t) = Cst_mt*max(0, Umt(t)-Umt(t-1))  + Cst_fc*max(0, Ufc(t)-Ufc(t-1));
    end
end
%% 目标函数每小时运行费用（目标函数1）
cost=[];
for t=1:24
    cost(t)=  Cch4*( (Pmt(t)/(L_gas*eta_mt(t))) + (Pfc(t)/(L_gas*eta_fc(t))) ) ...       % 燃料成本
                     +  Cm_mt*Pmt(t) + Cm_fc*Pfc(t)  +  Cm_pv*Ppv(t)...    % 维护成本
                     + Cm_wt*Pwt(t) + Cst(t)+(Crb(t)+1.5)*Px(t)+Cm_Eb*Pb(t);               
end
figure
plot(0:23,cost)
hold on
plot(0:23,costp)
title('每小时运行成本');
legend('本文调度方法','常规调度方法');
xlabel('时间/小时')
ylabel('费用/元')





