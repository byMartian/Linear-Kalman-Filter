function X_hat = KalmanFilter_FirstOrder(data,Q,R,x0,P0)
%KalmanFilter_FirstOrder һ�����Կ������˲���

N = length(data);

KalmanGain = zeros(N,1);
X_hat = zeros(N,1);
X_hat_pre = zeros(N,1);
P_covariance = zeros(N,1);
P_covariance_pre = zeros(N,1);

X_hat_pre(1) = x0;
P_covariance(1) = P0;

%ϵͳ����
T=0.05;%����ʱ����
A = 1;%״̬ת�ƾ���
H = 1; %��������
I= eye(1);

for k = 1:N-1
   %%%%%%%%%%%%%%%%%%% �������˲����� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_hat_pre(k+1) = A * X_hat(k);                                          %1)Ԥ�⣺����״̬����ֵ���������鹫ʽ�õ��Ĺ���ֵ�� ��3.36)
    P_covariance_pre(k+1) = A*P_covariance(k)*A' + Q;           %2)Ԥ�⣺�����飩Ԥ��Э������� ��3.39��
    
    KalmanGain(k+1) = P_covariance_pre(k+1)* H' * inv(H*P_covariance_pre(k+1)* H'+R);    %3)У��������Kalman������� ��3.38��
    %�������˲������Ĺ�ʽ����ǰ�Ĺ���ֵ(�˲����)=��һ�εĹ���ֵ+����������ϵ�� * (��ǰ�Ĳ���ֵ - ��һ�εĹ���ֵ)
    X_hat(k+1) = X_hat_pre(k+1) + KalmanGain(k+1) * (data(k+1) - H* X_hat_pre(k+1) );             %4)У��������״̬���� ���� �� ��3.37��
    P_covariance(k+1) = ( I-KalmanGain(k+1)*H ) * P_covariance_pre(k+1) ;                         %5)У�������Э������� ���£� ��3.40)
    
end

