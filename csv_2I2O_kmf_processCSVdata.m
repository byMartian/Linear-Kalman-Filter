%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%功能说明： 1）一阶线性卡尔曼滤波器+二阶线性卡尔曼滤波器的实例，用于处理[速度,加速度]状态估计；
%作者：chenzf2013@163.com
%参考：[1] DRCAN，https://www.bilibili.com/video/BV1hC4y1b7K7/?spm_id_from=333.788 
%      [2] 《黄小平, 王岩. 卡尔曼滤波原理及应用MATLAB仿真[M]. 电子工业出版社, 2015-有书签-扫描版》 式(3.36)到式(3.40)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 清除工作区的变量
clc;    
clear
%% (1.1)数据输入 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_raw = readmatrix('.\brake_thr.csv');%

% [kalmf,L,P] = kalman(sys,Q,R,N)
%% (1.2)数据提取
T=0.05; %数据采集周期
N=length(data_raw); %总的采样次数
Time_series = 0.05:0.05:(0.05*N);
 
 
speed_gps(:,1)=data_raw(:,2)/3.6;               %当前组合惯导的gps速度;m/s 
speed_current_log(:,1) = data_raw(:,10);      %当前踏歌ccu中所用滤波器记录的实际速度；
acc_imu(:,1) = data_raw(:,7)*9.8;              %惯导纵向加速度


%%  (2) 二阶线性卡尔曼滤波-2I2O-处理[gps车; imu纵向加速度]%%%%%%%%%%%%%%%%%%%%%%%%%
%% (2.1) 2I2O-参数: 噪声方差矩阵——可调参数
% 噪声；可调参数
% Q = [0.002,  0.00037;
%     -0.00037,  0.01]; %当前最好的参数
Q = [0.002,  0.00037;
    -0.00037,  0.01];
%过程噪声的协方差矩阵：两个状态值x=[x1 ; x2]，两个状态值的过程噪声分别为w=[w1; w2];
%过程噪声符合高斯分布p(w1)~(0,q11), p(w2)~(0,q22); 两个状态值的过程噪声组成的协方差矩阵为Q,  p(w)~(0,Q)；Q=[q11,q12; q21,q22 ]
% R = [0.01, 0.0; 
%         0.0,  1.8]; %当前最好的参数
R = [0.01, 0.0; 
    0.0,  1.4]; %当前最好的参数
%测量噪声的协方差矩阵；两个传感器的测量值z=[z1 ; z2]，两个传感器的测量噪声分别为v=[v1; v2];
%过程噪声符合高斯分布p(v1)~(0,r11), p(v2)~(0,r22); 两个状态值的过程噪声组成的协方差矩阵为R,  p(v)~(0,R)；R=[r11,r12; r21,r22 ]

%系统矩阵
T=0.05;%采样时间间隔
A = [1,T;   0,1];%状态转移矩阵
%B=[0.5*T*T;  T] ;%控制矩阵
H = [1, 0;  0,1]; %测量矩阵；
I = eye(2);

%% (2.2) 2I2O-初始化各矩阵
%定义各矩阵%%%%%%%%%%%%%%%%%%%
Z = cell(1,N); %测量值；
%预测
X_hat_pre = cell(1,N);%预测：先验估计值，依据先验公式得到的估计值；
P_covariance_pre = cell(1,N);%预测：先验误差协方差值；
%校正
X_hat = cell(1,N);%校正：后验估计值（卡尔曼滤波结果）；
KalmanGain = cell(1,N);%校正：卡尔曼增益；
P_covariance = cell(1,N);%校正：更新误差协方差；
%U = cell(1,N);%控制量

%初始化：各矩阵的值%%%%%%%%%%%%%%%%%%%
Z0 = [speed_gps(1,1); acc_imu(1,1)];%初始位移和速度
Z{1} = Z0;
X_hat{1} = Z0;%第一步的后验估计值（卡尔曼滤波结果）初始化：使用测量结果；
P_covariance{1} = [1, 0;     0, 1]; %第一步误差协方差矩阵化：单位矩阵来初始化；

%U{1} = 0;
P_covariance_pre{1} = P_covariance{1};
KalmanGain{1} = [0, 0;     0, 0]; 

%% (2.3)二阶线性卡尔曼滤波器-2I2O

for k = 1:N-1
      %%%%%%%%%%%%%%%%%%% 加速度传感器测量结果限幅 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Z{k+1} = [speed_gps(k+1,1); acc_imu(k+1,1)]; %迭代下一步：新的传感器测量值—[速度;加速度]
    %加速度增长限幅 0.15*g
    limit_acc_delta = 0.15 *9.8;
    limit_acc_delta_nav =  -0.15 *9.8;
    acc_delta = Z{k+1}(2,1) - X_hat{k}(2,1);
   if acc_delta > limit_acc_delta
        Z{k+1}(2,1) = limit_acc_delta + X_hat{k}(2,1);%测量的 加速度修正
    elseif acc_delta < limit_acc_delta_nav
        Z{k+1}(2,1) = limit_acc_delta_nav + X_hat{k}(2,1);%测量的 加速度修正
    else
       fprintf('###log### acc_delta[ %f]: %f\n', k,acc_delta);
    end
        
    %%%%%%%%%%%%%%%%%%% 卡尔曼滤波过程 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_hat_pre{k+1} = A * X_hat{k};                                          %1)预测：先验状态估计值，依据先验公式得到的估计值； （3.36)
    P_covariance_pre{k+1} = A*P_covariance{k}*A' + Q;           %2)预测：（先验）预测协方差矩阵 （3.39）
    
    KalmanGain{k+1} = P_covariance_pre{k+1}* H' * inv(H*P_covariance_pre{k+1}* H'+R);    %3)校正：计算Kalman增益矩阵； （3.38）
    %卡尔曼滤波器核心公式：当前的估计值(滤波结果)=上一次的估计值+卡尔曼增益系数 * (当前的测量值 - 上一次的估计值)
    X_hat{k+1} = X_hat_pre{k+1} + KalmanGain{k+1} * (Z{k+1} - H* X_hat_pre{k+1} );             %4)校正：后验状态估计 更新 ； （3.37）
    P_covariance{k+1} = ( I-KalmanGain{k+1}*H ) * P_covariance_pre{k+1} ;                         %5)校正：误差协方差矩阵 更新； （3.40)
    
end

%由cell变成矩阵，方便画图
XX_hat =cell2mat (X_hat);
ZZ =cell2mat (Z);

%% (2.3)+移动平均滤波
speed_2orderKF_mps= (XX_hat(1, :))';
acc_2orderKF = (XX_hat(2, :))';
speed_2orderKF_smooth = smooth(smooth(XX_hat(1, :),12) ,6);%两次：移动平均滤波+移动平均滤波; m/s
acc_2orderKF_smooth = smooth(XX_hat(2, :),10);%移动平均滤波;m/s2

acc_for_speed_difference = zeros(N,1);%由速度直接差分得到加速度；
acc_for_speed_difference(1)=acc_2orderKF_smooth(1);
for k= 1:N-1
    acc_for_speed_difference(k+1)= ( speed_2orderKF_smooth(k+1)-speed_2orderKF_smooth(k) )/0.05;
end
acc_for_speed_difference=smooth(acc_for_speed_difference,10);%结果移动平均滤波
    

%% (2.4)一阶线性卡尔曼滤波-SISO-处理gps车速

%%%%%%CCU中的参数%%%%%%%%%%%%%%%%%%%%%%%%
%  kalman_filter_spd_ = std::make_shared<math::KalmanFilter<>>(0.001,0.1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 噪声；可调参数
Q_oneDim_speed = 0.03; %当前最好的参数;
% Q_oneDim_speed = 0.001; %调节参数
R_oneDim_speedp =0.35;%当前最好的参数
% Q_oneDim_speed =0.1;
P0 = 1; %第一步误差协方差矩阵化：单位矩阵来初始化；
x0 = speed_gps(1) ;
speed_1orderKF= KalmanFilter_FirstOrder(speed_gps(:),Q_oneDim_speed,R_oneDim_speedp,x0,P0);
speed_1orderKF_kmh= 3.6*speed_1orderKF(:);

%%  (2.5)画图-速度(SISO+2I2O)  
figure('Name','速度滤波比较');
box on
plot(Time_series,ZZ(1,:)*3.6,'k');%测量速度
hold on
plot(Time_series, XX_hat(1, :)*3.6, '-r.',...
    'LineWidth',1,...
    'MarkerSize',12,...
    'MarkerFaceColor',[0.5,0.5,0.5]);% 二阶线性Kalman滤波器
hold on
plot(Time_series, speed_2orderKF_smooth(:)*3.6, '-c.',...
    'LineWidth',2,...
    'MarkerSize',12);% 二阶线性Kalman滤波器+移动平均滤波
hold on
plot(Time_series, speed_1orderKF_kmh, '-b.');% 一阶线性Kalman滤波器
hold on
plot(Time_series, speed_current_log(:), '-g.');% 当前ccu中的滤波速度
legend('gps测量速度', '二阶线性Kalman滤波器', '二阶线性Kalman滤波器+移动平均滤波','一阶线性Kalman滤波器','当前ccu中的速度-KF')
xlabel('时间 t(s)')
ylabel('速度 v(km/h)')
title('速度滤波比较');
hold off

%% (2.6)画图-加速度
figure('Name','加速度比较');
box on
plot(Time_series,ZZ(2,:),'k');%测量加速度
hold on
plot(Time_series, XX_hat(2, :), '-r.',...
    'LineWidth',1,...
    'MarkerSize',12,...
    'MarkerFaceColor',[0.5,0.5,0.5]);% 二阶Kalman 滤波
hold on
plot(Time_series,acc_2orderKF_smooth(:),'-c.',...
     'LineWidth',2,...
    'MarkerSize',12);%二阶线性Kalman滤波器+移动平均滤波
hold on
plot(Time_series,acc_for_speed_difference(:),'-g.');%由速度直接差分得到加速度
hold on

legend('imu测量加速度', '二阶线性Kalman滤波器', '二阶线性Kalman滤波器+移动平均滤波','由速度直接差分得到加速度')
xlabel('时间 t(s)')
ylabel('加速度 a(m/s2')
title('加速度比较');
hold off


%% (3)一阶线性卡尔曼滤波-SISO-处理imu俯仰角 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 噪声；可调参数
% Q_oneDim_pitch = 0.01; %当前最好的参数;
Q_oneDim_pitch = 0.01; %调节参数

% R_oneDim_pitch =0.05;%当前最好的参数
R_oneDim_pitch =0.05;

P0 = 1; %第一步误差协方差矩阵化：单位矩阵来初始化；
x0 = imu_pitch_angle(1) ;

imu_pitch_angle_kf(:)= KalmanFilter_FirstOrder(imu_pitch_angle(:),Q_oneDim_pitch,R_oneDim_pitch,x0,P0);
%移动平均滤波
imu_pitch_angle_kf_smooth = smooth(imu_pitch_angle_kf(:),12);%移动平均滤波

%% 画图-俯仰角相关
figure('Name','俯仰角比较');
box on
plot(Time_series,imu_pitch_angle(:,1) ,'k');%惯导俯仰角；度
hold on
plot(Time_series, imu_pitch_angle_filter(:,1) , '-r.',...
    'LineWidth',0.5,...
    'MarkerSize',6);% 惯导俯仰角滤波；度
hold on
plot(Time_series,imu_pitch_angle_kf(:),'-c.');%惯导俯仰角,KF滤波
hold on
plot(Time_series, imu_pitch_angle_kf_smooth(:) , '-g.',...
    'LineWidth',0.8,...
    'MarkerSize',8);% 惯导俯仰角滤波；度

legend('惯导俯仰角', '惯导俯仰角滤波','惯导俯仰角,KF滤波','惯导俯仰角,KF滤波+移动平均滤波')
xlabel('时间 t(s)')
ylabel('俯仰角 theta(度)')
title('惯导俯仰角比较');
hold off


%% (4)一阶线性卡尔曼滤波-SISO-处理坡度%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 噪声；可调参数
% Q_oneDim_slop = 0.03; %当前最好的参数;
Q_oneDim_slop = 0.01; %调节参数

% R_oneDim_slop =0.35;%当前最好的参数
R_oneDim_slop =0.05;

P0 = 1; %第一步误差协方差矩阵化：单位矩阵来初始化；
x0 = vehicle_slope_original(1) ;

vehicle_slope_1orderKF(:)= KalmanFilter_FirstOrder(vehicle_slope_original(:),Q_oneDim_slop,R_oneDim_slop,x0,P0);
vehicle_slope_1orderKF_smooth = smooth(vehicle_slope_1orderKF(:),12);%移动平均滤波

road_slope_for_imu_pitch = tan( deg2rad( imu_pitch_angle_kf_smooth(:)) );
road_slope_for_imu_pitch_smooth =smooth(road_slope_for_imu_pitch(:),6);
road_slope_for_imu_pitch_smooth_modify = road_slope_for_imu_pitch_smooth -0.031;
%% 画图-坡度相关
figure('Name','坡度比较');
box on
plot(Time_series,vehicle_slope_original(:,1) ,'k'); 
hold on
plot(Time_series, vehicle_slope_filter(:,1) , '-r.',...
    'LineWidth',0.5,...
    'MarkerSize',6); 
hold on
plot(Time_series,vehicle_slope_1orderKF(:),'-c.',...
     'LineWidth',0.5,...
    'MarkerSize',6); 
hold on
plot(Time_series, vehicle_slope_1orderKF_smooth(:), '-b.',...
        'LineWidth',1.0,...
    'MarkerSize',6);% 
hold on
plot(Time_series, road_slope_for_imu_pitch_smooth_modify(:), '-g.',...
        'LineWidth',1.5,...
    'MarkerSize',6);% 
legend('车辆坡度原始', 'CCU当前车辆坡度滤波', '车辆坡度KF滤波','车辆坡度KF滤波+移动平均滤波','来自惯导俯仰角换算得到+移动平均滤波+偏移修正')
xlabel('时间 t(s)')
ylabel('车辆坡度tan(theta)')
title('车辆坡度比较');
hold off


