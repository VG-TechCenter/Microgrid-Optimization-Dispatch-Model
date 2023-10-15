% Dung Beetle Optimizer: (DBO) (demo)
% Programmed by Jian-kai Xue    
% Updated 28 Nov.2022.                     
function [fMin, bestX, Convergence_curve] = IDBO(pop, M, c, d, dim, fobj)

P_percent = 0.2; % 执行滚球行为的蜣螂比例   
pNum = round(pop * P_percent); % 执行滚球行为的蜣螂数量   
lb = c.*ones(1,dim);    % 取值下限向量
ub = d.*ones(1,dim);    % 取值上限向量

% 种群初始化
% ★★改进1：Chebyshev映射★★
x = Chebyshev(pop, dim, ub, lb);
for i = 1 : pop
    fit(i) = fobj(x(i,:)) ;                       
end
pFit = fit;                       
pX = x; 
XX = pX;    
[fMin, bestI] = min(fit); % fMin是全局最优适应度值
bestX = x(bestI,:);       % bestX是与fMin对应的全局最优解

% 开始迭代，进行个体更新
for t = 1 : M 

[fmax, B] = max(fit);
worse = x(B,:);   
r2 = rand(1);

% 滚球行为蜣螂的位置更新（分为无障碍物模式和有障碍物模式）
for i = 1 : pNum    
    if(r2<0.9)
        % ①无障碍物模式
%         a = rand(1,1);
%         if (a>0.1)
%             a = 1;
%         else
%             a = -1;
%         end
%         x(i,:) = pX(i,:)+0.3*abs(pX(i,:)-worse)+a*0.1*(XX(i,:));
        % ★★改进2：黄金正弦策略★★
        R1 = rand()*2*pi;
        R2 = rand()*pi;
        tao = (sqrt(5)-1)/2;
        x1 = -pi+(1-tao)*2*pi;
        x2 = -pi+tao*2*pi; 
        x(i,:) = pX(i,:)*abs(sin(R1))+R2*sin(R1)*abs(x1*pX(i,:)-x2*x(i,:));
    else
        % ②有障碍物模式
        aaa = randperm(180,1);
        if (aaa==0 ||aaa==90 ||aaa==180)
            x(i,:) = pX(i,:);   
        end
        theta = aaa*pi/180;   
        x(i,:) = pX(i,:)+tan(theta).*abs(pX(i,:)-XX(i,:)); 
    end
    % 越界校正
    x(i,:) = Bounds(x(i,:), lb, ub);    
    fit(i) = fobj(x(i,:));
end

[fMMin, bestII] = min(fit);   % fMin是当前的最优适应度值
bestXX = x(bestII,:);         % bestXX是当前的最优解
R = 1-t/M;                         

Xnew1 = bestXX.*(1-R); 
Xnew2 = bestXX.*(1+R);                    
Xnew1 = Bounds(Xnew1, lb, ub);
Xnew2 = Bounds(Xnew2, lb, ub);

% 繁殖行为蜣螂的位置更新
for i = (pNum + 1) : 12               
    x(i,:) = bestXX+((rand(1,dim)).*(pX(i,:)-Xnew1)+(rand(1,dim)).*(pX(i,:)-Xnew2));
    x(i,:) = Bounds(x(i,:), Xnew1, Xnew2);
    fit(i) = fobj(x(i,:)) ;
end

Xnew11 = bestX.*(1-R); 
Xnew22 = bestX.*(1+R);                    
Xnew11 = Bounds(Xnew11, lb, ub);
Xnew22 = Bounds(Xnew22, lb, ub);

% 觅食行为小蜣螂的位置更新
for i = 13 : 19        
    x(i,:) = pX(i,:)+((randn(1)).*(pX(i,:)-Xnew11)+((rand(1,dim)).*(pX(i,:)-Xnew22)));
    x(i,:) = Bounds(x(i,:), lb, ub);
    fit(i) = fobj(x(i,:));
end

% 偷窃行为蜣螂的位置更新
% ★★改进3：位置更新动态权重系数★★
k1 = 1-((t)^3)/((M)^3);
k2 = ((t)^3)/((M)^3);
for j = 20 : pop       
    x(j,:) = k1*bestX+k2*randn(1,dim).*((abs((pX(j,:)-bestXX)))+(abs((pX(j,:)-bestX))))./2;
    x(j,:) = Bounds(x(j,:), lb, ub);
    fit(j) = fobj(x(j,:)) ;
end

% 更新个体最优值和全局最优值
XX = pX;
for i = 1 : pop 
    if (fit(i) < pFit(i))
        pFit(i) = fit(i);
        pX(i,:) = x(i,:);
    end

    if( pFit(i) < fMin)
        fMin = pFit(i);
        bestX = pX(i,:);
    end
end

% 迭代曲线记录
Convergence_curve(t) = fMin;
end

% 越界校正函数
function s = Bounds(s, Lb, Ub)
temp = s;
I = temp < Lb;
temp(I) = Lb(I);

J = temp > Ub;
temp(J) = Ub(J);
s = temp;

