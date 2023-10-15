function fun = objective(x)
%% 参数设置
parameter;
Pmt = x(1:24);            % 燃气轮机功率
Pfc = x(25:48);           % 燃料电池功率
Px = x(49:72);            % 可卸负荷                 
Pb = x(73:96);            % 蓄电池功率
% 预分配
eta_mt = zeros(1,24);%用来表示24个时段内燃气轮机的效率
eta_fc = zeros(1,24);%燃料电池24个时段的效率
%Pmth = zeros(1,24);%矩阵Pmth来表示燃气轮机热功率 
Umt = zeros(1,24);%定义矩阵Umt表示燃气轮机的启停状态
Ufc = zeros(1,24);%同理定义Ufc表示燃料电池的启停状态。

fun = 0;
g = [];
h = [];%约束方程组
%% 运行模型
for t=1:24
    % 燃气轮机热功率
    %下面这个是微型燃气轮机效率计算公式
    eta_mt(t) = 0.0753*(Pmt(t)/65)^3 - 0.3095*(Pmt(t)/65)^2 + 0.4174*(Pmt(t)/65) + 0.1068;
%     Pmth(t) = ((Pmt(t)*(1-eta_mt(t)-eta_l))/eta_mt(t))*eta_h*Coph;
     % 燃料电池功率
    eta_fc(t) = -0.0023*Pfc(t) + 0.674;
end
%% 电能平衡约束
%峰10-14   19-22  
%平7-10     14-19
%谷23-7
for t=1:24
    
      if Pb(t)<=0    % 对于微电网来说Pb(t)<0 代表充电,大于0代表放电，p1(t)代表负荷
        h = [h, Pmt(t)+Pfc(t)+Pwt(t)+Ppv(t)+Px(t)+Pb(t)/eta_bch-Pl(t)];%约束方程组
        %各种电力来源的总发电功率 + 电池提供的功率 - 负载需求功率 = 0
      elseif Pb(t)>0
        h = [h, Pmt(t)+Pfc(t)+Pwt(t)+Ppv(t)+Px(t)+Pb(t)*eta_bdis-Pl(t)];
        %之所以在[]中还要写h,是因为每加入一个新约束,需要把之前的约束方程也保留,才能形成约束方程组。
      end
    
end
 
%% 电池储能初始和最终状态相等约束
% % sum(x),假如x是一个矩阵，则sum(x)表示以矩阵x的每一列为对象，对一列数字进行求和
h = [h, sum(Pb)];
% 电储能约束
for t=1:24
    g=[g, Eb_init*(1-tau_b)-sum(Pb(1:t))-Eb_max];%最大储能量约束

% sum（pb(1:t)）表示访问了把pb从1到t个元素加起来
    g=[g, -(Eb_init*(1-tau_b)-sum(Pb(1:t))-Eb_min)];%最小储能量约束
end

%% 爬坡约束
for t=2:24
    % 燃气轮机爬坡约束
    g=[g, Pmt(t)-Pmt(t-1)-Rmt_up];
    g=[g, -(Pmt(t)-Pmt(t-1)-Rmt_down)];
    % 燃料电池爬坡约束
    g=[g, Pfc(t)-Pfc(t-1)-Rfc_up];
    g=[g, -(Pfc(t)-Pfc(t-1)-Rfc_down)];
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
%% 目标函数
cost=[];
for t=1:24
    cost(t)=  Cch4*( (Pmt(t)/(L_gas*eta_mt(t))) + (Pfc(t)/(L_gas*eta_fc(t))) ) ...       % 燃料成本
                     +  Cm_mt*Pmt(t) + Cm_fc*Pfc(t)  +  Cm_pv*Ppv(t)...    % 维护成本
                     + Cm_wt*Pwt(t) + Cst(t)+(Crb(t)+1.5)*Px(t)+Cm_Eb*Pb(t);               
end
fun=sum(cost);
%**********************罚函数处理***********************%
Big=10000;
small=0.01;
N=length(g);% 表示数列的长度，即有多少个元素
M=length(h);
G=0;
for n=1:N
    G=G+max(0, g(n))^2;
end
H=0;
for m=1:M
    H=H+max(  0, abs(h(m))-small  )^2;
end
%*************加入罚函数后的最终目标函数*************%
fun = fun+Big*(H+G);

end