%%
% 
% * R08323020 陳霖
% * R08323021 陸馨
% 


%% Q1
beta = 0.95;
delta = 0.5;
theta = 0.8;
alpha = 0.5;

A = 1;
k_star_old = ((1/beta+delta-1)/(alpha*A))^(1/(alpha-1));
c_star_old = A*k_star_old^alpha - delta*k_star_old;
k_star_old;
c_star_old;
%% capital before the shock
disp(k_star_old);
%% consumption before the shock
disp(c_star_old);

A = 1.1;
k_star_new = ((1/beta+delta-1)/(alpha*A))^(1/(alpha-1));
c_star_new = A*k_star_new^alpha - delta*k_star_new;
%% The steady states after the shock
%% capital before the shock
disp(k_star_new);
%% consumption before the shock
disp(c_star_new);



%% Q2, Q3
k_star = k_star_new;
c_star = c_star_new;

f2 = alpha*(alpha-1)*A*k_star^(alpha-2);
A_matrix = [1/beta -1; f2*c_star/theta 1-beta*f2*c_star/theta];

%% The linear approximate A
disp(A_matrix)

%%
e = eig(A_matrix);
[V,D] = eig(A_matrix);
%% eigenvalues
disp(e);
%% eigenvectors
disp(V);
%% Q4
k0_hat = k_star_old - k_star_new;
c0_hat = V(2)/V(1)*k0_hat;

k_vector = [k0_hat];
c_vector = [c0_hat];
kt_hat = k0_hat;
ct_hat = c0_hat;
for i = 6 : 30
    kt_hat = kt_hat * D(1);
    ct_hat = ct_hat * D(1);
    k_vector = [k_vector kt_hat];
    c_vector = [c_vector ct_hat];
end

k_vector = k_vector + k_star_new;
c_vector = c_vector + c_star_new;
k_vector = [k_star_old k_star_old k_star_old k_star_old k_star_old k_vector];
c_vector = [c_star_old c_star_old c_star_old c_star_old c_star_old c_vector];
output_vector = [];
investment_vector = [];

for i=1:5
   A = 1;
   output_vector(i) = A * k_vector(i)^alpha;
   investment_vector(i) = output_vector(i) - c_vector(i);    
end

for i=6:31
   A = 1.1;
   output_vector(i) = A * k_vector(i)^alpha;
   investment_vector(i) = output_vector(i) - c_vector(i);
end

y = 0:30;

%% Time Series of Capital
plot(y,k_vector,'-*b'); %線性，顏色，標記
set(gca,'XTick',[0:1:30]) %x軸範圍1-6，間隔1
xlabel('time')  %x軸座標描述
ylabel('k') %y軸座標描述
title('Time Series of Capital')
saveas(gcf,'k.png')

%% Time Series of Consumption
plot(y,c_vector,'-*b'); %線性，顏色，標記
set(gca,'XTick',[0:1:30]) %x軸範圍1-6，間隔1
xlabel('time')  %x軸座標描述
ylabel('c') %y軸座標描述
title('Time Series of Consumption')
saveas(gcf,'c.png')

%% Time Series of Output
plot(y,output_vector,'-*b'); %線性，顏色，標記
set(gca,'XTick',[0:1:30]) %x軸範圍1-6，間隔1
xlabel('time')  %x軸座標描述
ylabel('c') %y軸座標描述
title('Time Series of Output')
saveas(gcf,'output.png')

%% Time Series of Investment
plot(y,investment_vector,'-*b'); %線性，顏色，標記
set(gca,'XTick',[0:1:30]) %x軸範圍1-6，間隔1
xlabel('time')  %x軸座標描述
ylabel('c') %y軸座標描述
title('Time Series of Investment')
saveas(gcf,'investment.png')

%% Q5
theta = 0.1;

A = 1;
k_star_old = ((1/beta+delta-1)/(alpha*A))^(1/(alpha-1));
c_star_old = A*k_star_old^alpha - delta*k_star_old;
k_star_old;
c_star_old;

A = 1.1;
k_star_new = ((1/beta+delta-1)/(alpha*A))^(1/(alpha-1));
c_star_new = A*k_star_new^alpha - delta*k_star_new;
k_star_new;
c_star_new;

%
k_star = k_star_new;
c_star = c_star_new;

f2 = alpha*(alpha-1)*A*k_star^(alpha-2);
A_matrix = [1/beta -1; f2*c_star/theta 1-beta*f2*c_star/theta];
e = eig(A_matrix);
[V,D] = eig(A_matrix);
%% eigenvalues
disp(e);
%% eigenvectors
disp(V);

%%
k0_hat = k_star_old - k_star_new;
c0_hat = V(2)/V(1)*k0_hat;

k_vector = [k0_hat];
c_vector = [c0_hat];
kt_hat = k0_hat;
ct_hat = c0_hat;
for i = 6 : 30
    kt_hat = kt_hat * D(1);
    ct_hat = ct_hat * D(1);
    k_vector = [k_vector kt_hat];
    c_vector = [c_vector ct_hat];
end

k_vector = k_vector + k_star_new;
c_vector = c_vector + c_star_new;
k_vector = [k_star_old k_star_old k_star_old k_star_old k_star_old k_vector];
c_vector = [c_star_old c_star_old c_star_old c_star_old c_star_old c_vector];
output_vector = [];
investment_vector = [];

for i=1:5
   A = 1;
   output_vector(i) = A * k_vector(i)^alpha;
   investment_vector(i) = output_vector(i) - c_vector(i);    
end
for i=6:31
   A = 1.1;
   output_vector(i) = A * k_vector(i)^alpha;
   investment_vector(i) = output_vector(i) - c_vector(i);
end

y = 0:30;

%% Time Series of Capital (theta = 0.1)
plot(y,k_vector,'-*b'); %線性，顏色，標記
set(gca,'XTick',[0:1:30]) %x軸範圍1-6，間隔1
xlabel('time')  %x軸座標描述
ylabel('k') %y軸座標描述
title('Time Series of Capital (theta = 0.1)')
saveas(gcf,'k_5.png')

%% Time Series of Consumption (theta = 0.1)
plot(y,c_vector,'-*b'); %線性，顏色，標記
set(gca,'XTick',[0:1:30]) %x軸範圍1-6，間隔1
xlabel('time')  %x軸座標描述
ylabel('c') %y軸座標描述
title('Time Series of Consumption (theta = 0.1)')
saveas(gcf,'c_5.png')

%% Time Series of Output (theta = 0.1)
plot(y,output_vector,'-*b'); %線性，顏色，標記
set(gca,'XTick',[0:1:30]) %x軸範圍1-6，間隔1
xlabel('time')  %x軸座標描述
ylabel('c') %y軸座標描述
title('Time Series of Output (theta = 0.1)')
saveas(gcf,'output_5.png')


%% Time Series of Investment (theta = 0.1)
plot(y,investment_vector,'-*b'); %線性，顏色，標記
set(gca,'XTick',[0:1:30]) %x軸範圍1-6，間隔1
xlabel('time')  %x軸座標描述
ylabel('c') %y軸座標描述
title('Time Series of Investment (theta = 0.1)')
saveas(gcf,'investment_5.png')