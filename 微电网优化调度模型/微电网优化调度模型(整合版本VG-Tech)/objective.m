function fun = objective(x)
%% ��������
parameter;
Pmt = x(1:24);            % ȼ���ֻ�����
Pfc = x(25:48);           % ȼ�ϵ�ع���
Px = x(49:72);            % ��ж����                 
Pb = x(73:96);            % ���ع���
% Ԥ����
eta_mt = zeros(1,24);%������ʾ24��ʱ����ȼ���ֻ���Ч��
eta_fc = zeros(1,24);%ȼ�ϵ��24��ʱ�ε�Ч��
%Pmth = zeros(1,24);%����Pmth����ʾȼ���ֻ��ȹ��� 
Umt = zeros(1,24);%�������Umt��ʾȼ���ֻ�����ͣ״̬
Ufc = zeros(1,24);%ͬ����Ufc��ʾȼ�ϵ�ص���ͣ״̬��

fun = 0;
g = [];
h = [];%Լ��������
%% ����ģ��
for t=1:24
    % ȼ���ֻ��ȹ���
    %���������΢��ȼ���ֻ�Ч�ʼ��㹫ʽ
    eta_mt(t) = 0.0753*(Pmt(t)/65)^3 - 0.3095*(Pmt(t)/65)^2 + 0.4174*(Pmt(t)/65) + 0.1068;
%     Pmth(t) = ((Pmt(t)*(1-eta_mt(t)-eta_l))/eta_mt(t))*eta_h*Coph;
     % ȼ�ϵ�ع���
    eta_fc(t) = -0.0023*Pfc(t) + 0.674;
end
%% ����ƽ��Լ��
%��10-14   19-22  
%ƽ7-10     14-19
%��23-7
for t=1:24
    
      if Pb(t)<=0    % ����΢������˵Pb(t)<0 ������,����0����ŵ磬p1(t)������
        h = [h, Pmt(t)+Pfc(t)+Pwt(t)+Ppv(t)+Px(t)+Pb(t)/eta_bch-Pl(t)];%Լ��������
        %���ֵ�����Դ���ܷ��繦�� + ����ṩ�Ĺ��� - ���������� = 0
      elseif Pb(t)>0
        h = [h, Pmt(t)+Pfc(t)+Pwt(t)+Ppv(t)+Px(t)+Pb(t)*eta_bdis-Pl(t)];
        %֮������[]�л�Ҫдh,����Ϊÿ����һ����Լ��,��Ҫ��֮ǰ��Լ������Ҳ����,�����γ�Լ�������顣
      end
    
end
 
%% ��ش��ܳ�ʼ������״̬���Լ��
% % sum(x),����x��һ��������sum(x)��ʾ�Ծ���x��ÿһ��Ϊ���󣬶�һ�����ֽ������
h = [h, sum(Pb)];
% �索��Լ��
for t=1:24
    g=[g, Eb_init*(1-tau_b)-sum(Pb(1:t))-Eb_max];%�������Լ��

% sum��pb(1:t)����ʾ�����˰�pb��1��t��Ԫ�ؼ�����
    g=[g, -(Eb_init*(1-tau_b)-sum(Pb(1:t))-Eb_min)];%��С������Լ��
end

%% ����Լ��
for t=2:24
    % ȼ���ֻ�����Լ��
    g=[g, Pmt(t)-Pmt(t-1)-Rmt_up];
    g=[g, -(Pmt(t)-Pmt(t-1)-Rmt_down)];
    % ȼ�ϵ������Լ��
    g=[g, Pfc(t)-Pfc(t-1)-Rfc_up];
    g=[g, -(Pfc(t)-Pfc(t-1)-Rfc_down)];
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
%% Ŀ�꺯��
cost=[];
for t=1:24
    cost(t)=  Cch4*( (Pmt(t)/(L_gas*eta_mt(t))) + (Pfc(t)/(L_gas*eta_fc(t))) ) ...       % ȼ�ϳɱ�
                     +  Cm_mt*Pmt(t) + Cm_fc*Pfc(t)  +  Cm_pv*Ppv(t)...    % ά���ɱ�
                     + Cm_wt*Pwt(t) + Cst(t)+(Crb(t)+1.5)*Px(t)+Cm_Eb*Pb(t);               
end
fun=sum(cost);
%**********************����������***********************%
Big=10000;
small=0.01;
N=length(g);% ��ʾ���еĳ��ȣ����ж��ٸ�Ԫ��
M=length(h);
G=0;
for n=1:N
    G=G+max(0, g(n))^2;
end
H=0;
for m=1:M
    H=H+max(  0, abs(h(m))-small  )^2;
end
%*************���뷣�����������Ŀ�꺯��*************%
fun = fun+Big*(H+G);

end