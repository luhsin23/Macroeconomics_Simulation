%%
% 
% * R08323020 陳霖
% * R08323021 陸馨
% 

%% Q1
clear
A = 1.1;
alpha = 0.5;
delta = 0.5;
theta = 0.1;
beta = 0.95;
syms k_t
eqn = A * k_t^alpha + (1 - delta) * k_t == k_t;
k_bar = solve(eqn);
k_bar = eval(k_bar(2));
K = 0 : 0.002 : k_bar;      % 1 * 2421
n = length(K);              % 2421
%% The state variable
disp(K)

%% Q2, Q3
Y = A * K.^ alpha + (1 - delta) * K;
YY = ones(n, 1) * Y;        % 1 * 2421
KK = ones(n, 1) * K;        % 2421 * 2421
C = YY - KK';               % 2421 * 2421

% Check budget constraint: C>0
C_new = C;
C_new(C_new < 0) = NaN;
U = (C_new.^(1 - theta) - 1) / (1 - theta);

% Set the initial value function to all zeros
V = zeros(n ,1);            % 2421 * 1
VV = V * ones(1, n);        % 2421 * 2421
W = U + beta * VV;
% maximize over each k value
V = max(W)';                % 2421 * 1

%% Q4
epsilon = 1;
iter = 1;
while (epsilon > 10^ (-4)) && iter < 1000
VV = V * ones(1, n); 
W = U + beta * VV; 
% Check!
[TV, policy] = max(W);
V1 = TV';
epsilon = max(abs(V1 - V));
V = V1;
iter = iter + 1; 
end
%% Q4.2 policy rule
policy_rule = policy';
disp(policy_rule)
%% Q4.1 plot the value function
plot(V,'-.'); 
set(gca,'XTick', [0:200:n]) 
xlabel('k')
ylabel('v(k)')
title('value function')
saveas(gcf,'value_function.png')

%% Q5 plot dynamic path using value function iteration
A = 1;
k_star_old = ((1 / beta + delta - 1)/(alpha * A))^(1 / (alpha - 1));
% find the value in K less than but closest to k_star_old
index5 = floor(k_star_old/0.002) + 1;       % 410
k_index = [];
k_value = [];
index = index5;
for i = 6 : 31
    k_index_temp = policy(index);           % 444
    k_index = [k_index k_index_temp];
    k_value_temp = K(k_index_temp);         % 0.8860
    k_value = [k_value k_value_temp];
    index = k_index_temp; 
end
k_vector = ones(1, 6) * k_star_old;
k_vector = [k_vector k_value];

output_vector = [];
investment_vector = [];
A = 1;
c_star_old = A * k_star_old^alpha - delta * k_star_old;
output_old = A * k_star_old^alpha;
investment_old = output_old - c_star_old;
c_vector = ones(1, 5) * c_star_old;
output_vector = ones(1, 5) * output_old;
investment_vector = ones(1, 5) * investment_old;
for i = 6 : 31
   A = 1.1;
   output_vector(i) = A * k_vector(i)^alpha;
   c_vector(i) = output_vector(i) + (1 - delta) * k_vector(i) - k_vector(i+1);
   investment_vector(i) = output_vector(i) - c_vector(i);
end

% k, o correct
% check i, c
% a = [1 2 6 4 5];
% b = [1 3 4 5 6];
% plot(1:5, a);
% hold on
% plot(1:5,b);
% hold off


%% Q5.1 capital from value function iteration
plot(0 : 30, k_vector(1:31),'-*b');
set(gca,'XTick',[0 : 1 : 30]) 
xlabel('time')
ylabel('k')
title('capital')
saveas(gcf,'k_vfi.png')

%% Q5.2 consumption from value function iteration
plot(0 : 30, c_vector,'-*b');
set(gca,'XTick',[0 : 1 : 30]) 
xlabel('time')
ylabel('c')
title('consumption')
saveas(gcf,'c_vfi.png')

%% Q5.3 output from value function iteration
plot(0 : 30, output_vector,'-*b');
set(gca,'XTick',[0 : 1 : 30]) 
xlabel('time')
ylabel('output')
title('output')
saveas(gcf,'o_vfi.png')

%% Q5.4 investment from value function iteration
plot(0 : 30, investment_vector,'-*b');
set(gca,'XTick',[0 : 1 : 30]) 
xlabel('time')
ylabel('investment')
title('investment')
saveas(gcf,'i_vfi.png')


%% Q6 comparasion of value function iteration and linear approximation
m_linear = readtable('M.csv');
k_linear = m_linear.Var2;
c_linear = m_linear.Var3;
output_linear = m_linear.Var4;
investment_linear = m_linear.Var5;


%% Q6.1 capital comparasion
a1 = plot(0 : 30, k_vector(1:31),'-*b', 'Color', 'blue');
hold on
a2 = plot(0 : 30, k_linear,'-*b', 'Color', 'red');
hold off
set(gca,'XTick',[0 : 1 : 30]) 
xlabel('time')
ylabel('capital')
title('capital')
M1 = 'value function iteraton';
M2 = 'linear approximation';
legend([a1;a2], M1, M2,'Location','best');
saveas(gcf,'k_compare.png')

%% Q6.2 consumption comparasion
a1 = plot(0 : 30, c_vector,'-*b', 'Color', 'blue');
hold on
a2 = plot(0 : 30, c_linear,'-*b', 'Color', 'red');
hold off
set(gca,'XTick',[0 : 1 : 30]) 
xlabel('time')
ylabel('consumption')
title('consumption')
M1 = 'value function iteraton';
M2 = 'linear approximation';
legend([a1;a2], M1, M2,'Location','best');
saveas(gcf,'c_compare.png')

%% Q6.3 output comparasion
a1 = plot(0 : 30, output_vector,'-*b', 'Color', 'blue');
hold on
a2 = plot(0 : 30, output_linear,'-*b', 'Color', 'red');
hold off
set(gca,'XTick',[0 : 1 : 30]) 
xlabel('time')
ylabel('output')
title('output')
M1 = 'value function iteraton';
M2 = 'linear approximation';
legend([a1;a2], M1, M2,'Location','best');
saveas(gcf,'o_compare.png')

%% Q6.4 investment comparasion
a1 = plot(0 : 30, investment_vector,'-*b', 'Color', 'blue');
hold on
a2 = plot(0 : 30, investment_linear,'-*b', 'Color', 'red');
hold off
set(gca,'XTick',[0 : 1 : 30]) 
xlabel('time')
ylabel('investment')
title('investment')
M1 = 'value function iteraton';
M2 = 'linear approximation';
legend([a1;a2], M1, M2,'Location','best');
saveas(gcf,'i_compare.png')