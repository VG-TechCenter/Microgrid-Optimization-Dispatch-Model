clear;
close all;
clc;

SearchAgents_no = 30; % ��Ⱥ����
Function_name = 'F2'; % ���Ժ��������� F1 �� F23 
Max_iteration = 500;  % ����������
cnt_max = 30;         % ���д���
[lb, ub, dim, fobj] = Get_Functions_details(Function_name); % ������ѡ���Ժ�������ϸ��Ϣ
Curve_DBO = zeros(1, Max_iteration);
Curve_IDBO = zeros(1, Max_iteration);

for cnt = 1 : cnt_max
    [DBO_Best_score(cnt), DBO_Best_pos(cnt,:), DBO_curve] = DBO(SearchAgents_no, Max_iteration, lb, ub, dim, fobj);      % ԭʼDBO�㷨
    [IDBO_Best_score(cnt), IDBO_Best_pos(cnt,:), IDBO_curve] = IDBO(SearchAgents_no, Max_iteration, lb, ub, dim, fobj);  % �Ľ����DBO�㷨
    Curve_DBO = Curve_DBO + DBO_curve;
    Curve_IDBO = Curve_IDBO + IDBO_curve;
end

% ��������
Curve_DBO = Curve_DBO/cnt_max;
Curve_IDBO = Curve_IDBO/cnt_max;

% ��׼��
std_DBO = std(DBO_Best_score);
std_IDBO = std(IDBO_Best_score);

% ���ֵ
worst_DBO = max(DBO_Best_score);
worst_IDBO = max(IDBO_Best_score);

% ����ֵ
best_DBO = min(DBO_Best_score);
best_IDBO = min(IDBO_Best_score);

% ����ֵ
mean_DBO = mean(DBO_Best_score);
mean_IDBO = mean(IDBO_Best_score);

% ��ͼ
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
title([num2str(Function_name), '��Ѱ������']);
xlabel('Iteration');
ylabel('Best score obtained so far');
grid on
box on
legend('DBO','IDBO')

disp(['������', num2str(Function_name)]);
disp(['DBO�����ֵ: ', num2str(worst_DBO), ', ����ֵ: ', num2str(best_DBO), ', ƽ��ֵ: ', num2str(mean_DBO), ', ��׼��: ', num2str(std_DBO)]);
disp(['IDBO�����ֵ: ', num2str(worst_IDBO), ', ����ֵ: ', num2str(best_IDBO), ', ƽ��ֵ: ', num2str(mean_IDBO), ', ��׼��: ', num2str(std_IDBO)]);




