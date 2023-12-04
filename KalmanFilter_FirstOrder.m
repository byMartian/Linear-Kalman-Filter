function X_hat = KalmanFilter_FirstOrder(data,Q,R,x0,P0)
%KalmanFilter_FirstOrder 一阶线性卡尔曼滤波器

N = length(data);

KalmanGain = zeros(N,1);
X_hat = zeros(N,1);
X_hat_pre = zeros(N,1);
P_covariance = zeros(N,1);
P_covariance_pre = zeros(N,1);

X_hat_pre(1) = x0;
P_covariance(1) = P0;

%系统矩阵
T=0.05;%采样时间间隔
A = 1;%状态转移矩阵
H = 1; %测量矩阵
I= eye(1);

for k = 1:N-1
   %%%%%%%%%%%%%%%%%%% 卡尔曼滤波过程 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_hat_pre(k+1) = A * X_hat(k);                                          %1)预测：先验状态估计值，依据先验公式得到的估计值； （3.36)
    P_covariance_pre(k+1) = A*P_covariance(k)*A' + Q;           %2)预测：（先验）预测协方差矩阵 （3.39）
    
    KalmanGain(k+1) = P_covariance_pre(k+1)* H' * inv(H*P_covariance_pre(k+1)* H'+R);    %3)校正：计算Kalman增益矩阵； （3.38）
    %卡尔曼滤波器核心公式：当前的估计值(滤波结果)=上一次的估计值+卡尔曼增益系数 * (当前的测量值 - 上一次的估计值)
    X_hat(k+1) = X_hat_pre(k+1) + KalmanGain(k+1) * (data(k+1) - H* X_hat_pre(k+1) );             %4)校正：后验状态估计 更新 ； （3.37）
    P_covariance(k+1) = ( I-KalmanGain(k+1)*H ) * P_covariance_pre(k+1) ;                         %5)校正：误差协方差矩阵 更新； （3.40)
    
end

