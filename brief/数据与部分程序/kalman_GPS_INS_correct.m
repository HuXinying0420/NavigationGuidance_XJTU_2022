function [ E_att, E_pos, E_vn, Xk, Pk ] = kalman_GPS_INS_correct( Xk_1, Qk_1, Pk_1, F, G, H, ts ,Zk,Rk )
%kalman�˲���λ�����

%����˵��
%����
%Zk   ����λ����������� ��Ϊ�˲�������
%Xk_1   k-1ʱ��״̬����������ֵ��
%Qk_1   ϵͳ�����������
%Rk     ���������������
%Pk_1   ���ƾ���������
%F,G,H  ����ϵͳ��״̬ת�ƾ�����������Ͳ�������
%ts     ��������
%���
%E_pos     λ���������Ԥ��ֵ
%E-vn      �ٶ��������Ԥ��ֵ
%E-att     ��̬�������Ԥ��ֵ
%Xk        kʱ��״̬����ֵ����
%Pk        kʱ�̹��ƾ���������

%���������ʼ��
E_att = zeros(3,1);
E_pos = E_att;
E_vn = E_att;

%����ϵͳ��ɢ��
%���в�������HΪ������
Phikk_1 = eye(18,18)+F*ts;
Gammak_1 = (eye(18,18)+ts*F/2)*G*ts;
Hk=H;

if nargin<9    %������״̬���ƣ�����GPS��������Ƶ�ʱ�INS��
    Xk = Phikk_1*Xk_1;
    Pk = Phikk_1*Pk_1*Phikk_1'+Gammak_1*Qk_1*Gammak_1';
else            % �в���ʱ�˲�
     %�������˲�
    Pkk_1 = Phikk_1*Pk_1*Phikk_1'+Gammak_1*Qk_1*Gammak_1';     %�������Ԥ�����
    Kk = Pkk_1*Hk'/(Hk*Pkk_1*Hk'+Rk);    %�������

    %Pk��ʽ����֤Pk�ĶԳ���
    temp=eye(18)-Kk*Hk;  %��ʱ����
    Pk=temp*Pkk_1*temp'+Kk*Rk*Kk';
    
    Xkk_1= Phikk_1*Xk_1;         %״̬һ��Ԥ��
    Xk = Xkk_1+Kk*(Zk-Hk*Xkk_1); %״̬����ֵ
 
    
    E_att= Xk(1:3);   %�������
    E_vn= Xk(4:6);
    E_pos= Xk(7:9);
end

end

