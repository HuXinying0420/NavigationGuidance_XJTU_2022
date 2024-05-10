function [ E_att, E_pos, E_vn, Xk, Pk ] = kalman_GPS_INS_correct( Xk_1, Qk_1, Pk_1, F, G, H, ts ,Zk,Rk )
%kalman滤波，位置组合

%参数说明
%输入
%Zk   量测位置误差向量， 作为滤波器输入
%Xk_1   k-1时刻状态向量（估计值）
%Qk_1   系统噪声方差矩阵
%Rk     测量噪声方差矩阵
%Pk_1   估计均方误差矩阵
%F,G,H  连续系统的状态转移矩阵，噪声矩阵和测量矩阵
%ts     迭代步长
%输出
%E_pos     位置误差向量预测值
%E-vn      速度误差向量预测值
%E-att     姿态误差向量预测值
%Xk        k时刻状态估计值向量
%Pk        k时刻估计均方误差矩阵

%输出参数初始化
E_att = zeros(3,1);
E_pos = E_att;
E_vn = E_att;

%连续系统离散化
%其中测量矩阵H为常矩阵
Phikk_1 = eye(18,18)+F*ts;
Gammak_1 = (eye(18,18)+ts*F/2)*G*ts;
Hk=H;

if nargin<9    %仅进行状态递推，由于GPS采样更新频率比INS低
    Xk = Phikk_1*Xk_1;
    Pk = Phikk_1*Pk_1*Phikk_1'+Gammak_1*Qk_1*Gammak_1';
else            % 有测量时滤波
     %卡尔曼滤波
    Pkk_1 = Phikk_1*Pk_1*Phikk_1'+Gammak_1*Qk_1*Gammak_1';     %均方误差预测矩阵
    Kk = Pkk_1*Hk'/(Hk*Pkk_1*Hk'+Rk);    %增益矩阵

    %Pk公式，保证Pk的对称性
    temp=eye(18)-Kk*Hk;  %临时变量
    Pk=temp*Pkk_1*temp'+Kk*Rk*Kk';
    
    Xkk_1= Phikk_1*Xk_1;         %状态一步预测
    Xk = Xkk_1+Kk*(Zk-Hk*Xkk_1); %状态估计值
 
    
    E_att= Xk(1:3);   %输出储存
    E_vn= Xk(4:6);
    E_pos= Xk(7:9);
end

end

