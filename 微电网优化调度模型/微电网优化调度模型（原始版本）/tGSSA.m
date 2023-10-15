%_________________________________________________________________________%
%% [1]��ΰ��,����.����Ӧt�ֲ���ƽ����ҸĽ�����ȸ�����㷨����Ӧ��[J/OL].΢����ѧ������:1-8[2021-12-17].https://doi.org/10.19304/J.ISSN1000-7180.2020-0026.
%_________________________________________________________________________%
function [Best_pos,Best_score,curve]=tGSSA(fobj,dim,lb,ub,Max_iter,pop)

ST = 0.6;%Ԥ��ֵ
PD = 0.7;%�����ߵı��У�ʣ�µ��Ǽ�����
SD = 0.2;%��ʶ����Σ����ȸ�ı���

PDNumber = round(pop*PD); %����������
SDNumber = round(SD*pop);%��ʶ����Σ����ȸ����
if(max(size(ub)) == 1)
   ub = ub.*ones(1,dim);
   lb = lb.*ones(1,dim);  
end

%��Ⱥ��ʼ��
X0=initialization(pop,dim,ub,lb);
X = X0;
%�����ʼ��Ӧ��ֵ
fitness = zeros(1,pop);
for i = 1:pop
   fitness(i) =  fobj(X(i,:));
end
 [fitness, index]= sort(fitness);%����
BestF = fitness(1);
WorstF = fitness(end);
GBestF = fitness(1);%ȫ��������Ӧ��ֵ
for i = 1:pop
    X(i,:) = X0(index(i),:);
end
curve=zeros(1,Max_iter);
GBestX = X(1,:);%ȫ������λ��
X_new = X;
for i = 1: Max_iter
    
    BestF = fitness(1);
    WorstF = fitness(end);

    
    R2 = rand(1);
   for j = 1:PDNumber
      if(R2<ST)
          %% �Ľ��㣺�ƽ����ҸĽ�
          %�ƽ����Ҳ���
          a = 2 - 2*i/Max_iter;
          b = 1;
          r=rand();
          r1=(2*pi)*r;
          r2=r*pi;
          gold=double((sqrt(5)-1)/2);      % �ƽ�ָ���, around 0.618
          x1=a+(1-gold)*(b-a);            %�ƽ�ָ�ϵ��x1
          x2=a+gold*(b-a);                %�ƽ�ָ�ϵ��x2
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
          %����-1��1�������
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
   %�߽����
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
   %����λ��
   for j=1:pop
    fitness_new(j) = fobj(X_new(j,:));
   end
   %% �Ľ��㣺 ����Ӧt�ֲ�
   w1 = 0.5;w2 =0.1;
   p = w1-w2*(Max_iter-i)/Max_iter;%��̬ѡ�����
   for j = 1:pop
       if(p<rand)
          Temp = X_new(j,:) +X_new(j,:)*trnd(i); %���ڵ���������t�ֲ�����
           %�߽紦��
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
    %�������
   [fitness, index]= sort(fitness);%����
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



