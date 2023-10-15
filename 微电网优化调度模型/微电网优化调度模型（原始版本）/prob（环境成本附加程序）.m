%得到多目标问题的解（目标函数2）
function [y,c] = prob(x)  %c=1则x为非可行解

[c,y(1)] = fitness(x);  %得出粒子的适应度  %%运行成本赋值给y（1）
%% 目标函数2:环境保护成本
C_DE_en=0;C_grid_en=0;C_MT_en=0;
for i=1:144
    if i>72&&i<97
      C_DE_en=C_DE_en+(0.023*680+6*0.306+8*10.09)*x(i); 
    elseif  i>96&&i<121
     C_MT_en=C_MT_en+((0.023*889+6*1.8+8*1.6)*6*x(i));
    end  
end
for i=121:144
    if x(i)>0
        C_grid_en=C_grid_en+(0.023*724+6*0.0036+8*0.2)*x(i);
    else
        C_grid_en=C_grid_en+0;
    end
end
C_en= C_DE_en+ C_MT_en+ C_grid_en;
y(2) =C_en; %环境保护成本赋值给y（2）

end