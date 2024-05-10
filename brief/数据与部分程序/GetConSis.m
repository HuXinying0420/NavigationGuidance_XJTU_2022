function [ F,G,H ] = GetConSis( vn, pos, quat, Fn, Tg, Ta )
%�õ�����ϵͳ��״̬����ϵ������F,G,H
%dx = F*x + G*w
%z = H*x + v
%����˵��
%pos    ins���λ����������Ϊ�˲���ϵͳ����
%vn     ins����ٶ���������Ϊ�˲���ϵͳ����
%Fn     ins�Ľ������룬����ϵ�±�����������Ϊ�˲�������
%quat   ins�����̬��Ԫ������Ϊ�˲�������
%Ta     ���ٶȼ����Ư�����ʱ��
%Tg     ���������Ư�����ʱ��
%���
%F      ����ϵͳ״̬ת�ƾ���18*18
%G      ����ϵͳ��������18*9
%H      ����ϵͳ��������3*18

glvs;

%��������ʼ��
Re = glv.Re;   %���򳤰뾶
f = glv.f;     %�������
wie= glv.wie;  %������ת���ٶ�

%�������ٶ�
Ve=vn(1);Vn=vn(2);Vu=vn(3);
%����λ��
L=pos(1);h=pos(3);


fe = Fn(1);fn = Fn(2);fu = Fn(3);
Rm = Re*(1-2*f+3*f*sin(L)^2);
Rn = Re*(1+f*sin(L)^2);

q = quat;
Cnb = [1-2*(q(3)^2+q(4)^2),     2*(q(2)*q(3)-q(1)*q(4)), 2*(q(2)*q(4)+q(1)*q(3));
    2*(q(2)*q(3)+q(1)*q(4)), 1-2*(q(2)^2+q(4)^2),     2*(q(3)*q(4)-q(1)*q(2));
    2*(q(2)*q(4)-q(1)*q(3)), 2*(q(3)*q(4)+q(1)*q(2)), 1-2*(q(2)^2+q(3)^2)];

%����ϵͳ״̬ת���� F ��ʱ�����
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
%����ϵͳ����������
G = zeros(18,9);
G(1:3,1:3)   = Cnb;
G(13:15,4:6) = eye(3,3);
G(16:18,7:9) = eye(3,3);
%����ϵͳ���������
H = zeros(3,18);
H(1,7)= 1;
H(2,8)= 1;
H(3,9)= 1;

end

