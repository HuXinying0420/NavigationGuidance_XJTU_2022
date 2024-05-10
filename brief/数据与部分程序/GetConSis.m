function [ F,G,H ] = GetConSis( vn, pos, quat, Fn, Tg, Ta )
%得到连续系统的状态方程系数矩阵F,G,H
%dx = F*x + G*w
%z = H*x + v
%参数说明
%pos    ins输出位置向量，作为滤波器系统参数
%vn     ins输出速度向量，作为滤波器系统参数
%Fn     ins的解算输入，导航系下比力向量，作为滤波器参数
%quat   ins输出姿态四元数，作为滤波器参数
%Ta     加速度计误差漂移相关时间
%Tg     陀螺仪误差漂移相关时间
%输出
%F      连续系统状态转移矩阵18*18
%G      连续系统噪声矩阵18*9
%H      连续系统测量矩阵3*18

glvs;

%各参数初始化
Re = glv.Re;   %地球长半径
f = glv.f;     %地球扁率
wie= glv.wie;  %地球自转角速度

%东北天速度
Ve=vn(1);Vn=vn(2);Vu=vn(3);
%导航位置
L=pos(1);h=pos(3);


fe = Fn(1);fn = Fn(2);fu = Fn(3);
Rm = Re*(1-2*f+3*f*sin(L)^2);
Rn = Re*(1+f*sin(L)^2);

q = quat;
Cnb = [1-2*(q(3)^2+q(4)^2),     2*(q(2)*q(3)-q(1)*q(4)), 2*(q(2)*q(4)+q(1)*q(3));
    2*(q(2)*q(3)+q(1)*q(4)), 1-2*(q(2)^2+q(4)^2),     2*(q(3)*q(4)-q(1)*q(2));
    2*(q(2)*q(4)-q(1)*q(3)), 2*(q(3)*q(4)+q(1)*q(2)), 1-2*(q(2)^2+q(3)^2)];

%连续系统状态转换阵 F 的时间更新
F = zeros(18,18);
F(1,2)  = wie*sin(L)+Ve*tan(L)/(Rn+h);
F(1,3)  = -(wie*cos(L)+Ve/(Rn+h));
F(1,5)  = -1/(Rm+h);
F(1,9)  = Vn/(Rm+h)^2;
F(2,1)  = -(wie*sin(L)+Ve*tan(L)/(Rn+h));
F(2,3)  = -Vn/(Rm+h);
F(2,4)  = 1/(Rn+h);
F(2,7)  = -wie*sin(L);
F(2,9)  = -Ve/(Rn+h)^2;
F(3,1)  = wie*cos(L)+Ve/(Rn+h);
F(3,2)  = Vn/(Rm+h);
F(3,4)  = tan(L)/(Rn+h);
F(3,7)  = wie*cos(L)+Ve*(sec(L)^2)/(Rn+h);
F(3,9)  = -Ve*tan(L)/(Rn+h)^2;
F(4,2)  = -fu;
F(4,3)  = fn;
F(4,4)  = Vn*tan(L)/(Rm+h)-Vu/(Rm+h);
F(4,5)  = 2*wie*sin(L)+Ve*tan(L)/(Rn+h);
F(4,6)  = -(2*wie*cos(L)+Ve/(Rn+h));
F(4,7)  = 2*wie*cos(L)*Vn+Ve*Vn*sec(L)^2/(Rn+h)+2*wie*sin(L)*Vu;
F(4,9)  = (Ve*Vu-Ve*Vn*tan(L))/(Rn+h)^2;
F(5,1)  = fu;
F(5,3)  = -fe;
F(5,4)  = -2*(wie*sin(L)+Ve*tan(L)/(Rn+h));
F(5,5)  = -Vu/(Rm+h);
F(5,6)  = -Vn/(Rm+h);
F(5,7)  = -(2*wie*cos(L)+Ve*(sec(L)^2)/(Rn+h))*Ve;
F(5,9)  = (Ve^2*tan(L)+Vn*Vu)/(Rn+h)^2;
F(6,1)  = -fn;
F(6,2)  = fe;
F(6,4)  = 2*(wie*cos(L)+Ve/(Rn+h));
F(6,5)  = 2*Vn/(Rm+h);
F(6,7)  = -2*Ve*wie*sin(L);
F(6,9)  = -(Vn^2+Ve^2)/(Rn+h)^2;
F(7,5)  = 1/(Rm+h);
F(8,4)  = 1/((Rn+h)*cos(L));
F(8,7)  = Ve*tan(L)/((Rn+h)*cos(L));
F(9,6)  = 1;
F(1:3,10:12) = Cnb;
F(1:3,13:15) = Cnb;
F(4:6,16:18) = Cnb;
F(13,13) = -1/Tg(1);
F(14,14) = -1/Tg(2);
F(15,15) = -1/Tg(3);
F(16,16) = -1/Ta(1);
F(17,17) = -1/Ta(2);
F(18,18) = -1/Ta(3);
%连续系统输入矩阵更新
G = zeros(18,9);
G(1:3,1:3)   = Cnb;
G(13:15,4:6) = eye(3,3);
G(16:18,7:9) = eye(3,3);
%连续系统量测阵更新
H = zeros(3,18);
H(1,7)= 1;
H(2,8)= 1;
H(3,9)= 1;

end

