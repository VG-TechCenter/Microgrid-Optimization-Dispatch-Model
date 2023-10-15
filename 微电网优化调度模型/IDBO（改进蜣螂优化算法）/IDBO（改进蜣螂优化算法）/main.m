clear;
close all;
clc;

SearchAgents_no = 30; % 种群数量
Function_name = 'F2'; % 测试函数名：从 F1 到 F23 
Max_iteration = 500;  % 最大迭代次数
cnt_max = 30;         % 运行次数
[lb, ub, dim, fobj] = Get_Functions_details(Function_name); % 加载所选测试函数的详细信息
Curve_DBO = zeros(1, Max_iteration);
Curve_IDBO = zeros(1, Max_iteration);

for cnt = 1 : cnt_max
    [DBO_Best_score(cnt), DBO_Best_pos(cnt,:), DBO_curve] = DBO(SearchAgents_no, Max_iteration, lb, ub, dim, fobj);      % 原始DBO算法
    [IDBO_Best_score(cnt), IDBO_Best_pos(cnt,:), IDBO_curve] = IDBO(SearchAgents_no, Max_iteration, lb, ub, dim, fobj);  % 改进后的DBO算法
    Curve_DBO = Curve_DBO + DBO_curve;
    Curve_IDBO = Curve_IDBO + IDBO_curve;
end

% 迭代曲线
Curve_DBO = Curve_DBO/cnt_max;
Curve_IDBO = Curve_IDBO/cnt_max;

% 标准差
std_DBO = std(DBO_Best_score);
std_IDBO = std(IDBO_Best_score);

% 最差值
worst_DBO = max(DBO_Best_score);
worst_IDBO = max(IDBO_Best_score);

% 最优值
best_DBO = min(DBO_Best_score);
best_IDBO = min(IDBO_Best_score);

% 最优值
mean_DBO = mean(DBO_Best_score);
mean_IDBO = mean(IDBO_Best_score);

% 画图
figure(1);
func_plot(Function_name);
title('Test function')
xlabel('x_1');
ylabel('x_2');
zlabel([Function_name,'( x_1 , x_2 )'])
grid off

figure(2)
semilogy(Curve_DBO,'Color','b','linewidth', 1.5)
hold on
semilogy(Curve_IDBO,'Color','r','linewidth', 1.5)
title([num2str(Function_name), '的寻优曲线']);
xlabel('Iteration');
ylabel('Best score obtained so far');
grid on
box on
legend('DBO','IDBO')

disp(['函数：', num2str(Function_name)]);
disp(['DBO：最差值: ', num2str(worst_DBO), ', 最优值: ', num2str(best_DBO), ', 平均值: ', num2str(mean_DBO), ', 标准差: ', num2str(std_DBO)]);
disp(['IDBO：最差值: ', num2str(worst_IDBO), ', 最优值: ', num2str(best_IDBO), ', 平均值: ', num2str(mean_IDBO), ', 标准差: ', num2str(std_IDBO)]);




