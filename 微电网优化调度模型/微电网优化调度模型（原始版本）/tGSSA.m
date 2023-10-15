%_________________________________________________________________________%
%% [1]张伟康,刘升.自适应t分布与黄金正弦改进的麻雀搜索算法及其应用[J/OL].微电子学与计算机:1-8[2021-12-17].https://doi.org/10.19304/J.ISSN1000-7180.2020-0026.
%_________________________________________________________________________%
function [Best_pos,Best_score,curve]=tGSSA(fobj,dim,lb,ub,Max_iter,pop)

ST = 0.6;%预警值
PD = 0.7;%发现者的比列，剩下的是加入者
SD = 0.2;%意识到有危险麻雀的比重

PDNumber = round(pop*PD); %发现者数量
SDNumber = round(SD*pop);%意识到有危险麻雀数量
if(max(size(ub)) == 1)
   ub = ub.*ones(1,dim);
   lb = lb.*ones(1,dim);  
end

%种群初始化
X0=initialization(pop,dim,ub,lb);
X = X0;
%计算初始适应度值
fitness = zeros(1,pop);
for i = 1:pop
   fitness(i) =  fobj(X(i,:));
end
 [fitness, index]= sort(fitness);%排序
BestF = fitness(1);
WorstF = fitness(end);
GBestF = fitness(1);%全局最优适应度值
for i = 1:pop
    X(i,:) = X0(index(i),:);
end
curve=zeros(1,Max_iter);
GBestX = X(1,:);%全局最优位置
X_new = X;
for i = 1: Max_iter
    
    BestF = fitness(1);
    WorstF = fitness(end);

    
    R2 = rand(1);
   for j = 1:PDNumber
      if(R2<ST)
          %% 改进点：黄金正弦改进
          %黄金正弦参数
          a = 2 - 2*i/Max_iter;
          b = 1;
          r=rand();
          r1=(2*pi)*r;
          r2=r*pi;
          gold=double((sqrt(5)-1)/2);      % 黄金分割率, around 0.618
          x1=a+(1-gold)*(b-a);            %黄金分割系数x1
          x2=a+gold*(b-a);                %黄金分割系数x2
          X_new(j,:) = X(j,:)*abs(sin(r1)) - r2.*sin(r1).*abs(x1.*GBestX-x2*X(j,:)); 
      else
          X_new(j,:) = X(j,:) + randn().*ones(1,dim);
      end     
   end
   for j = PDNumber+1:pop
%        if(j>(pop/2))
        if(j>(pop - PDNumber)/2 + PDNumber)
          X_new(j,:)= randn(1,dim).*exp((X(end,:) - X(j,:))/j^2);
       else
          %产生-1，1的随机数
          A = ones(1,dim);
          for a = 1:dim
            if(rand()>0.5)
                A(a) = -1;
            end
          end 
          AA = A'*inv(A*A');     
          X_new(j,:)= X(1,:) + abs(X(j,:) - X(1,:)).*AA';
       end
   end
   Temp = randperm(pop);
   SDchooseIndex = Temp(1:SDNumber); 
   for j = 1:SDNumber
       if(fitness(SDchooseIndex(j))>BestF)
           X_new(SDchooseIndex(j),:) = X(1,:) + randn().*abs(X(SDchooseIndex(j),:) - X(1,:));
       elseif(fitness(SDchooseIndex(j))== BestF)
           K = 2*rand() -1;
           X_new(SDchooseIndex(j),:) = X(SDchooseIndex(j),:) + K.*(abs( X(SDchooseIndex(j),:) - X(end,:))./(fitness(SDchooseIndex(j)) - fitness(end) + 10^-8));
       end
   end
   %边界控制
   for j = 1:pop
       for a = 1: dim
           if(X_new(j,a)>ub)
               X_new(j,a) =ub(a);
           end
           if(X_new(j,a)<lb)
               X_new(j,a) =lb(a);
           end
       end
   end 
   %更新位置
   for j=1:pop
    fitness_new(j) = fobj(X_new(j,:));
   end
   %% 改进点： 自适应t分布
   w1 = 0.5;w2 =0.1;
   p = w1-w2*(Max_iter-i)/Max_iter;%动态选择概率
   for j = 1:pop
       if(p<rand)
          Temp = X_new(j,:) +X_new(j,:)*trnd(i); %基于迭代次数的t分布变异
           %边界处理
           Temp(Temp>ub) = ub(Temp>ub);
           Temp(Temp<lb) = lb(Temp<lb);
           fitvalue = fobj(Temp);
           if(fitvalue <fitness_new(j))
               X_new(j,:) = Temp;
               fitness_new(j) = fitvalue;
           end 
       end
   end
   
   
   for j = 1:pop
    if(fitness_new(j) < GBestF)
       GBestF = fitness_new(j);
        GBestX = X_new(j,:);   
    end
   end
   X = X_new;
   fitness = fitness_new;
    %排序更新
   [fitness, index]= sort(fitness);%排序
   BestF = fitness(1);
   WorstF = fitness(end);
   for j = 1:pop
      X(j,:) = X(index(j),:);
   end
   curve(i) = GBestF;
end
Best_pos =GBestX;
Best_score = curve(end);
end



