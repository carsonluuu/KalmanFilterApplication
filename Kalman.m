%% ========Part 1.1: Parameters setup========
clear;
tic
close all;
amp = 10; 
freq = 2*pi*0.1;
time = 30;
acc_f = 200;
gps_f = 5;
acc_t = 1/200;
gps_t = 1/5;
dt = 0.005; %1/acc_f
monte = 1;
%monte = 5000;
t_acc = 0:acc_t:time; %1x6001
t_gps = 0:gps_t:time; %1x151
t = 0:dt:time; %1x6001 used as acc_f as sample time

% means and variances are below sd -> sqrt(var)
mean_w = 0;
var_w = 0.0004;  sd_w = 0.02;
mean_ba = 0;
var_ba = 0.01;  sd_ba = 0.1;

mean_P0 = 0;
var_P0 = 100;  sd_P0 = 10;
mean_V0 = 100;
var_V0 = 1;  sd_V0 = 1; 
var_n1 = 1;  sd_n1 = 1; 
var_n2 = 0.04^2;  sd_n2 = 0.04;

%initialize some arrays to store 
n = zeros(2,length(t));%n stands for eta
%accelerometer model values
% a_c = zeros(1,length(t));
v_c = zeros(1,length(t));
p_c = zeros(1,length(t));
%truth values
a_t = zeros(1,length(t));
v_t = zeros(1,length(t));
p_t = zeros(1,length(t));
%% ========Part 1.2: Parameters Randomlize========
for j_monte = 1:1:monte
rng('shuffle');
v0 = 100 + sd_V0*randn;
p0 = 2;%sd_P0*randn;
n(1,:) = sd_n1*randn(1,length(t));
n(2,:) = sd_n2*randn(1,length(t));
ba = mean_ba + sd_ba*randn;
ba = ones(1,length(t))*ba;
w = mean_w + randn(1,length(t_acc))*sd_w;
%% ========Part 1.3: Model Matrix Assigned========
phi = [1 dt -0.5*dt*dt;0 1 -dt;0 0 1];
gama = -[0.5*dt*dt; dt; 0];
M0 = [var_V0 0 0; 0 var_P0 0; 0 0 var_ba];
W0 = var_w;
H = [1 0 0;0 1 0];
V0 = [var_n1 0;0 var_n2];
x_bar = [0.5;0.6;0.08];
%x_bar = [p0;v0;ba(1)]; 
%x_bar = [0.5;sd_V0*randn;sd_ba*randn];   %not sure...
%% ========Part 2.1: True Value Based on Motion========
a = amp*sin(freq*t);
v = v0 + (amp/freq) - (amp/freq)*cos(freq*t);
p = p0 + (v0 + (amp/freq))*t - (amp/freq^2)*sin(freq*t);
%% ========Part 2.2: Accelerometer Singal========
a_c = a + ba + w; %model equation
v_c(1) = mean_V0;
p_c(1)= mean_P0;

for k = 2:length(t_acc) 
    v_c(k) = v_c(k-1) + a_c(k-1)*acc_t;
    p_c(k) = p_c(k-1) + v_c(k-1)*acc_t + a_c(k-1)*0.5*acc_t^2;
end
%% ========Part 2.3: Sampled into Accelerometer for Truth======== 
for j = 1:length(t_acc) %store true value from motion
 a_t(j) = a(j);
 v_t(j) = v(j);
 p_t(j) = p(j);
end

%error 
Aerr_truth_acc = a_t - a_c;
Verr_truth_acc = v_t - v_c; 
Perr_truth_acc = p_t - p_c; 
%% Part 2.4: GPS Measurement
gps_p = p_t(1:acc_f/gps_f:length(t)) + n(1,1:acc_f/gps_f:length(t));
gps_v = v_t(1:acc_f/gps_f:length(t)) + n(2,1:acc_f/gps_f:length(t));
dp = p_t(1:acc_f/gps_f:length(t)) - p_c(1:acc_f/gps_f:length(t));
dv = v_t(1:acc_f/gps_f:length(t)) - v_c(1:acc_f/gps_f:length(t));

z(1,:) = dp + n(1,1:acc_f/gps_f:length(t));
z(2,:) = dv + n(2,1:acc_f/gps_f:length(t));
%% ========Part 3: Kalman Filter Implementation=======
kkk = 1;
M = M0;
for kk = 1:length(t)
    if (rem(kk-1, 40) == 0)
%             K = M0*H'/(H*M0*H'+V0);
%             r = z(:,kkk) - H*x_bar;
%             x_hat = x_bar + K*r;
%             P0 = (eye(3)-K*H)*M0*(eye(3)-K*H)' + K*V0*K';
% 
%             x_bar = phi*x_hat;
%                 xbar(:,kk) = x_bar; 
%             M = phi*P0*phi'+gama*W0*gama';
%             M0 = M;
%         
%                 xhat(:,kkk) = x_hat;

            P = inv(inv(M)+H'*inv(V0)*H);
                P_store(kkk) = {P};
            K = P*H'*inv(V0);
                K_store(kkk) = {K};
            x_hat = x_bar +K*(z(:,kkk)-H*x_bar);
                res(:,kkk) = z(:,kkk)-H*x_bar;
            x_bar = phi*x_hat;
            M = phi*P*phi'+gama*W0*gama';

                xhat(:,kkk) = x_hat;
                kkk = kkk + 1;
    else
        x_bar = phi*x_bar; 
        M = phi*M*phi'+gama*W0*gama';
    end
end

error(1,:) = xhat(1,:) -  (p_t(1:40:length(t)) - p_c(1:40:length(t)));
error(2,:) = xhat(2,:) +  v_c(1:40:length(t)) - v_t(1:40:length(t));
error(3,:) = xhat(3,:) - ba(1:40:length(t));

P_matrice = cell2mat(P_store(1:151));
K_matrice = cell2mat(K_store(1:151));
xhat_monte(j_monte) = {xhat};
res_monte(j_monte) = {res};
error_p(j_monte) = {error(1,:)};
error_v(j_monte) = {error(2,:)};
error_ba(j_monte) = {error(3,:)};
end

%% 
sum_p = 0;
sum_v = 0;
sum_ba = 0;
sum_P_ave = 0;
sum_orth = 0;
sum_res_orth = 0;
for num = 1:monte
    sum_p = sum_p + cell2mat(error_p(num));
    ave_p = sum_p/monte;
    sum_v = sum_v + cell2mat(error_v(num));
    ave_v = sum_v/monte;
    sum_ba = sum_ba + cell2mat(error_ba(num));
    ave_ba = sum_ba/monte;
end
Pave={zeros(3,3)};
for num_t = 1:151
    for num = 1:monte
        mat_p = cell2mat(error_p(num));
        mat_v = cell2mat(error_v(num));
        mat_ba = cell2mat(error_ba(num));
        e_l = [mat_p;mat_v;mat_ba]; 
        e_ave = [ave_p;ave_v;ave_ba];
        sum_P_ave = sum_P_ave + (e_l(:,num_t) - e_ave(:,num_t))*...
            (e_l(:,num_t) - e_ave(:,num_t))';
        P_ave = sum_P_ave./(monte-1);
    end
    Pave(num_t)={P_ave};
end

P_ave_mat = cell2mat(Pave(1:151));
P_store_mat = cell2mat(P_store(1:151));

for num_t = 1:151
    for num = 1:monte
        xhat_orth = cell2mat(xhat_monte(num));
        e_l = [mat_p;mat_v;mat_ba]; 
        e_ave = [ave_p;ave_v;ave_ba];
        sum_orth = sum_orth + (e_l(:,num_t) - e_ave(:,num_t))*...
            xhat_orth(:,num_t)';
        orth = sum_orth./(monte);
    end
    Orth(num_t)={orth};
end

Orth_mat = cell2mat(Orth(1:151));

for num = 1:monte
        res_orth = cell2mat(res_monte(num));
        sum_res_orth = sum_res_orth + res_orth(:,100)*...
            res_orth(:,50)';
        ensemble_orth = sum_res_orth./(monte);
end


%% 
% figure(1)
% subplot(311)
% plot(t,a,'g')
% ylabel('$Acceleration$ $(m/s^2)$','FontSize',12,'Interpreter','Latex');
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% title('\bf{Acceleration from Truth Model}','FontSize',11, 'Interpreter','Latex');
% subplot(312)
% plot(t,v,'g')
% ylabel('$Velocity$ $(m/s)$','FontSize',12,'Interpreter','Latex');
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% title('\bf{Velocity from Truth model}','FontSize',11, 'Interpreter','Latex');
% subplot(313)
% plot(t,p,'g')
% ylabel('$Position$ $(m)$','FontSize',12,'Interpreter','Latex');
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% title('\bf{Position from Truth Model}','FontSize',11, 'Interpreter','Latex');
% 
% figure(2)
% subplot(311)
% plot(t,Aerr_truth_acc)
% ylabel('$Acceleration$ $(m/s^2)$','FontSize',12,'Interpreter','Latex');
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% title('\bf{Acceleration Error between Accerometer and Truth Model}','FontSize',11, 'Interpreter','Latex');
% subplot(312)
% plot(t,Verr_truth_acc,'m')
% ylabel('$Velocity$ $(m/s)$','FontSize',12,'Interpreter','Latex');
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% title('\bf{Velocity Error between Accerometer and Truth Model}','FontSize',11, 'Interpreter','Latex');
% subplot(313)
% plot(t,Perr_truth_acc,'m')
% ylabel('$Position$ $(m)$','FontSize',12,'Interpreter','Latex');
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% title('\bf{Position Error between Accerometer and Truth Model}','FontSize',11, 'Interpreter','Latex');
% 
% figure(3)
% subplot(211)
% plot(t_gps,z(2,:),'r')
% ylabel('$Velocity$ $(m/s)$','FontSize',12,'Interpreter','Latex');
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% title('\bf{Velocity Error between Accerometer and Truth Model for Measurement}'...
%     ,'FontSize',10, 'Interpreter','Latex');
% subplot(212)
% plot(t_gps,z(1,:),'r')
% ylabel('$Position$ $(m)$','FontSize',12,'Interpreter','Latex');
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% title('\bf{Position Error between Accerometer and Truth Model for Measurement}'...
%     ,'FontSize',10, 'Interpreter','Latex');
% 
% figure(4)
% subplot(211)
% plot(t_gps,n(2,1:acc_f/gps_f:length(t)),'y')
% ylabel('$Velocity$ $(m/s)$','FontSize',12,'Interpreter','Latex');
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% title('\bf{Velocity Error between GPS output and truth model(noises in GPS measurements)}'...
%     ,'FontSize',8, 'Interpreter','Latex');
% subplot(212)
% plot(t_gps,n(1,1:acc_f/gps_f:length(t)),'y')
% ylabel('$Position$ $(m)$','FontSize',12,'Interpreter','Latex');
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% title('\bf{Position Error between GPS output and truth model(noises in GPS measurements)}'...
%     ,'FontSize',8, 'Interpreter','Latex');
% 
% figure(5)
% subplot(211)
% plot(t,Verr_truth_acc,'m')
% hold on
% plot(t_gps,z(2,:))
% ylabel('$Velocity$ $(m/s)$','FontSize',12,'Interpreter','Latex');
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% title('\bf{Velocity Error from Accerometer bias}'...
%     ,'FontSize',12, 'Interpreter','Latex');
% h51 = legend ( 'Accelerometer(model)','GPS(measurement)');
% set(h51,'FontSize',14);
% 
% subplot(212)
% plot(t,Perr_truth_acc,'m')
% hold on
% plot(t_gps,z(1,:))
% ylabel('$Position$ $(m)$','FontSize',12,'Interpreter','Latex');
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% title('\bf{Position Error from Accerometer bias}'...
%     ,'FontSize',12, 'Interpreter','Latex');
% h51 = legend ( 'Accelerometer(model)','GPS(measurement)');
% set(h51,'FontSize',14);
% 
% figure(61)
% plot(t,Perr_truth_acc,'m')
% hold on
% plot(t_gps,z(1,:))
% hold on
% plot(t_gps,xhat(1,:),'g-')
% ylabel('$Position$ $(m)$','FontSize',12,'Interpreter','Latex');
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% title('\bf{Position Error from Accerometer bias}'...
%     ,'FontSize',12, 'Interpreter','Latex');
% h61 = legend ( 'Accelerometer(model)','GPS(measurement)','Kalman');
% set(h61,'FontSize',14);
% 
% figure(62)
% plot(t,Verr_truth_acc,'m')
% hold on
% plot(t_gps,z(2,:))
% hold on
% plot(t_gps,xhat(2,:),'g-')
% ylabel('$Velocity$ $(m/s)$','FontSize',12,'Interpreter','Latex');
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% title('\bf{Velocity Error from Accerometer bias}'...
%     ,'FontSize',12, 'Interpreter','Latex');
% h62 = legend ( 'Accelerometer(model)','GPS(measurement)','Kalman');
% set(h62,'FontSize',14);
% 
% figure(7)
% plot(t,ba,'m.','Markersize',5)
% hold on
% plot(t_gps,xhat(3,:))
% ylabel('$Bias$ $(m/s^2)$','FontSize',12,'Interpreter','Latex');
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% title('\bf{Estimated Bias and Actual Bias}'...
%     ,'FontSize',12, 'Interpreter','Latex');
% h7 = legend ( 'Actual Bias','Estimated Bias');
% set(h7,'FontSize',14);

figure(81)
plot(t_gps,error(1,:),t_gps,sqrt(P_matrice(1,1:3:453)),t_gps,-sqrt(P_matrice(1,1:3:453)))
ylabel('$Position$ $(m)$','FontSize',12,'Interpreter','Latex');
xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
title('\bf{Position Error from Filtered and Modeled in Accerometer output and Truth Model}'...
    ,'FontSize',9, 'Interpreter','Latex');
ylim([-1.2 1.2])
h81 = legend ('Position Error','\sigma','-\sigma');
set(h81,'FontSize',14);

figure(82)
plot(t_gps,error(2,:),t_gps,sqrt(P_matrice(2,2:3:453)),t_gps,-sqrt(P_matrice(2,2:3:453)))
ylabel('$Velocity$ $(m/s)$','FontSize',12,'Interpreter','Latex');
xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
title('\bf{Velocity Error from Filtered and Modeled in Accerometer output and Truth Model}'...
    ,'FontSize',9, 'Interpreter','Latex');
ylim([-0.1 0.1])
h82 = legend ('Velocity Error','\sigma','-\sigma');
set(h82,'FontSize',14);

figure(83)
plot(t_gps,error(3,:),t_gps,sqrt(P_matrice(3,3:3:453)),t_gps,-sqrt(P_matrice(3,3:3:453)))
ylabel('$bias$ $(m)$','FontSize',12,'Interpreter','Latex');
xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
title('\bf{Bias Error between Filtered and Actual Value}'...
    ,'FontSize',9, 'Interpreter','Latex');
ylim([-0.12 0.12])
h83 = legend ('Bias Error','\sigma','-\sigma');
set(h83,'FontSize',14);

% figure(9)
% 
% subplot(321)
% plot(t_gps,(K_matrice(1,1:2:302)))
% y1 = ylabel('$K_{11}$','FontSize',10,'Interpreter','Latex','rot',0);
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% set(y1,'Units','Normalized','Position',[0.1,0.5,0]);
% 
% subplot(322)
% plot(t_gps,(K_matrice(1,2:2:302)))
% y2 = ylabel('$K_{12}$','FontSize',10,'Interpreter','Latex','rot',0);
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% set(y2,'Units','Normalized','Position',[0.1,0.5,0]);
% 
% subplot(323)
% plot(t_gps,(K_matrice(2,1:2:302)))
% y3 = ylabel('$K_{21}$','FontSize',10,'Interpreter','Latex','rot',0);
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% set(y3,'Units','Normalized','Position',[0.1,0.5,0]);
% 
% subplot(324)
% plot(t_gps,(K_matrice(2,2:2:302)))
% y4 = ylabel('$K_{22}$','FontSize',10,'Interpreter','Latex','rot',0);
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% set(y4,'Units','Normalized','Position',[0.1,0.5,0]);
% 
% subplot(325)
% plot(t_gps,K_matrice(3,1:2:302))
% y5 = ylabel('$K_{31}$','FontSize',10,'Interpreter','Latex','rot',0);
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% set(y5,'Units','Normalized','Position',[0.1,0.5,0]);
% 
% subplot(326)
% plot(t_gps,K_matrice(3,2:2:302))
% y6 = ylabel('$K_{32}$','FontSize',10,'Interpreter','Latex','rot',0);
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% set(y6,'Units','Normalized','Position',[0.1,0.5,0]);

% figure(10)
% 
% subplot(311)
% plot(t_gps,0.01*e_l(1,:))
% ylabel('$Position$ $(m/s^2)$','FontSize',12,'Interpreter','Latex');
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% title('\bf{Position Averaged error from Monte Carlo Simulation}','FontSize',11, 'Interpreter','Latex');
% ylim([-9e-3 9e-3])
% subplot(312)
% plot(t_gps,0.05*e_l(2,:))
% ylabel('$Velocity$ $(m/s)$','FontSize',12,'Interpreter','Latex');
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% title('\bf{Velocity Averaged error from Monte Carlo Simulation}','FontSize',11, 'Interpreter','Latex');
% ylim([-9e-3 9e-3])
% subplot(313)
% plot(t_gps,0.1*e_l(3,:))
% ylabel('$Bais$ $(m)$','FontSize',12,'Interpreter','Latex');
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% title('\bf{Bias Averaged error from Monte Carlo Simulation}','FontSize',11, 'Interpreter','Latex');
% ylim([-9e-3 9e-3])
% figure(11)
% 
% subplot(331)
% plot(t_gps,P_ave_mat(1,1:3:453) -P_store_mat(1,1:3:453))
% y1 = ylabel('$P_{11}$','FontSize',10,'Interpreter','Latex','rot',0);
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% set(y1,'Units','Normalized','Position',[0.1,0.5,0]);
% 
% subplot(332)
% plot(t_gps,P_ave_mat(1,2:3:453) -P_store_mat(1,2:3:453))
% y2 = ylabel('$P_{12}$','FontSize',10,'Interpreter','Latex','rot',0);
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% set(y2,'Units','Normalized','Position',[0.1,0.5,0]);
% 
% subplot(333)
% plot(t_gps,P_ave_mat(1,3:3:453) -P_store_mat(1,3:3:453)-0.004)
% y3 = ylabel('$P_{13}$','FontSize',10,'Interpreter','Latex','rot',0);
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% set(y3,'Units','Normalized','Position',[0.1,0.5,0]);
% 
% subplot(334)
% plot(t_gps,P_ave_mat(2,1:3:453) -P_store_mat(2,1:3:453))
% y4 = ylabel('$P_{21}$','FontSize',10,'Interpreter','Latex','rot',0);
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% set(y4,'Units','Normalized','Position',[0.1,0.5,0]);
% 
% subplot(335)
% plot(t_gps,P_ave_mat(2,2:3:453) -P_store_mat(2,2:3:453))
% y5 = ylabel('$P_{22}$','FontSize',10,'Interpreter','Latex','rot',0);
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% set(y5,'Units','Normalized','Position',[0.1,0.5,0]);
% 
% subplot(336)
% plot(t_gps,P_ave_mat(2,3:3:453) - P_store_mat(2,3:3:453))
% y6 = ylabel('$P_{23}$','FontSize',10,'Interpreter','Latex','rot',0);
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% set(y6,'Units','Normalized','Position',[0.1,0.5,0]);
% 
% subplot(337)
% plot(t_gps,P_ave_mat(3,1:3:453) -P_store_mat(3,1:3:453))
% y7 = ylabel('$P_{31}$','FontSize',10,'Interpreter','Latex','rot',0);
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% set(y7,'Units','Normalized','Position',[0.1,0.5,0]);
% 
% subplot(338)
% plot(t_gps,P_ave_mat(3,2:3:453) - P_store_mat(3,2:3:453))
% y8 = ylabel('$P_{32}$','FontSize',10,'Interpreter','Latex','rot',0);
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% set(y8,'Units','Normalized','Position',[0.1,0.5,0]);
% 
% subplot(339)
% plot(t_gps,P_ave_mat(3,3:3:453)-P_store_mat(3,3:3:453))
% y9 = ylabel('$P_{33}$','FontSize',10,'Interpreter','Latex','rot',0);
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% set(y9,'Units','Normalized','Position',[0.1,0.5,0]);

% figure(12)
% 
% subplot(331)
% plot(t_gps,Orth_mat(1,1:3:453)
% y1 = ylabel('$O_{11}$','FontSize',10,'Interpreter','Latex','rot',0);
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% set(y1,'Units','Normalized','Position',[0.1,0.5,0]);
% 
% subplot(332)
% plot(t_gps,Orth_mat(1,2:3:453))
% y2 = ylabel('$O_{12}$','FontSize',10,'Interpreter','Latex','rot',0);
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% set(y2,'Units','Normalized','Position',[0.1,0.5,0]);
% 
% subplot(333)
% plot(t_gps,Orth_mat(1,3:3:453))
% y3 = ylabel('$O_{13}$','FontSize',10,'Interpreter','Latex','rot',0);
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% set(y3,'Units','Normalized','Position',[0.1,0.5,0]);
% 
% subplot(334)
% plot(t_gps,Orth_mat(2,1:3:453))
% y4 = ylabel('$O_{21}$','FontSize',10,'Interpreter','Latex','rot',0);
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% set(y4,'Units','Normalized','Position',[0.1,0.5,0]);
% 
% subplot(335)
% plot(t_gps,Orth_mat(2,2:3:453))
% y5 = ylabel('$O_{22}$','FontSize',10,'Interpreter','Latex','rot',0);
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% set(y5,'Units','Normalized','Position',[0.1,0.5,0]);
% 
% subplot(336)
% plot(t_gps,Orth_mat(2,3:3:453))
% y6 = ylabel('$O_{23}$','FontSize',10,'Interpreter','Latex','rot',0);
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% set(y6,'Units','Normalized','Position',[0.1,0.5,0]);
% 
% subplot(337)
% plot(t_gps,Orth_mat(3,1:3:453))
% y7 = ylabel('$O_{31}$','FontSize',10,'Interpreter','Latex','rot',0);
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% set(y7,'Units','Normalized','Position',[0.1,0.5,0]);
% 
% subplot(338)
% plot(t_gps,Orth_mat(3,2:3:453))
% y8 = ylabel('$O_{32}$','FontSize',10,'Interpreter','Latex','rot',0);
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% set(y8,'Units','Normalized','Position',[0.1,0.5,0]);
% 
% subplot(339)
% plot(t_gps,Orth_mat(3,3:3:453)
% y9 = ylabel('$O_{33}$','FontSize',10,'Interpreter','Latex','rot',0);
% xlabel('$Time$ $(s)$','FontSize',10,'Interpreter','Latex');
% set(y9,'Units','Normalized','Position',[0.1,0.5,0]);
%toc