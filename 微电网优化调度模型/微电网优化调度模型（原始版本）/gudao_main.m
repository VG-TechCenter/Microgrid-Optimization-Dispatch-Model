clc;
clear;
close all;
global costp Ppv Pwt %costp���ܴ���ɱ�����,Ppv��Pwt��ʾ����ͷ������繦��
%% �㷨����
parameter;                %parameter�ؼ����������������Ĳ���,ʹ�������ƶԺ�����ɼ�,�����ں������ھͿ���ֱ��ʹ����Щ����������
nVar=4*24;                % Number of Decision Variables  ���ڱ�ʾ�Ż�����ľ��߱�������
VarMin=[ones(1,24)*Pmt_min, ones(1,24)*Pfc_min, ones(1,24)*Px_min, ones(1,24)*Pb_min];%�Ż��������½�����
VarMax=[ones(1,24)*Pmt_max, ones(1,24)*Pfc_max, ones(1,24)*Px_max, ones(1,24)*Pb_max];%�Ż��������Ͻ�����
MaxIt=100;                % Maximum Number of Iterations ����������
nPop=500;                 % Population Size (Swarm Size) ����Ⱥ���ģ

%% ����
[ bestPosition, fitValue ] =SSA(@objective,nVar,VarMin,VarMax,MaxIt,nPop);
x=bestPosition;
Pmt = x(1:24);            % ȼ���ֻ�����
Pfc = x(25:48);           % ȼ�ϵ�ع���
Px = x(49:72);            % ��ж����               
Pb = x(73:96);            % ���ع���
t=1:24;

%% ����������Ԥ��
figure
plot(t,Ppv,'-')
title('�����������');
xlabel('ʱ��/Сʱ')
ylabel('����/kw')
%% ��������������Ԥ��
figure
plot(t,Pwt,'-')
title('������������');
xlabel('ʱ��/Сʱ')
ylabel('����/kw')
%% �����ƽ����
figure
hold on 
Pb_po=max(Pb,0);
Pb_ne=min(Pb,0);
positive=[Pmt', Pfc', Pb_po',Px',Ppv',Pwt'];
negative=[ Pb_ne'];
bar(positive,'stack');
bar(negative,'stack');
plot(t, Pl, 'ok-');
title('��ƽ��');
legend('Pmtȼ��','Pfcȼ���','Pbdis��طŵ�','Px��ȥ����','Ppv�������','Pwt���' ,'Pbch��س��','Pl�ܸ���');
grid on
hold off
xlabel('ʱ��/Сʱ')
ylabel('����/kw')


figure
plot(t, Pl, 'ok-');
set(gca,'xtick',1:24)
title('�ܸ���');
grid on
xlabel('ʱ��/Сʱ')
ylabel('����/kw')
%% ��������ֳ������
figure
plot(t,Pmt,'ok-')
hold on
plot(t,Pfc,'-*')
hold on
plot(t,Px,'-')
hold on
plot(t,Pb,'-.')
legend('Pmtȼ������','Pfcȼ��س���','Px��ȥ����','Pb��س���');
title('����ͼ');
xlabel('ʱ��/Сʱ')
ylabel('����/kw')
%% ����ÿСʱ���з���
% Ԥ����
eta_mt = zeros(1,24);
eta_fc = zeros(1,24);
Pmth = zeros(1,24);
Umt = zeros(1,24);
Ufc = zeros(1,24);


%% ����ģ��
for t=1:24
    % ȼ���ֻ��ȹ���
    %���������΢��ȼ���ֻ�Ч�ʼ��㹫ʽ
    eta_mt(t) = 0.0753*(Pmt(t)/65)^3 - 0.3095*(Pmt(t)/65)^2 + 0.4174*(Pmt(t)/65) + 0.1068;
    Pmth(t) = ((Pmt(t)*(1-eta_mt(t)-eta_l))/eta_mt(t))*eta_h*Coph;
     % ȼ�ϵ�ع���
    eta_fc(t) = -0.0023*Pfc(t) + 0.674;
end



%% ��ͣ 
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
%% Ŀ�꺯��ÿСʱ���з��ã�Ŀ�꺯��1��
cost=[];
for t=1:24
    cost(t)=  Cch4*( (Pmt(t)/(L_gas*eta_mt(t))) + (Pfc(t)/(L_gas*eta_fc(t))) ) ...       % ȼ�ϳɱ�
                     +  Cm_mt*Pmt(t) + Cm_fc*Pfc(t)  +  Cm_pv*Ppv(t)...    % ά���ɱ�
                     + Cm_wt*Pwt(t) + Cst(t)+(Crb(t)+1.5)*Px(t)+Cm_Eb*Pb(t);               
end
figure
plot(0:23,cost)
hold on
plot(0:23,costp)
title('ÿСʱ���гɱ�');
legend('���ĵ��ȷ���','������ȷ���');
xlabel('ʱ��/Сʱ')
ylabel('����/Ԫ')





