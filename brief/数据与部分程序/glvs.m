%ȫ�ֱ���
global glv
glv.Re=6378160;%����뾶(������)
glv.f=1/298.3;%�������
glv.e=sqrt(2*glv.f-glv.f^2);%ƫ����
glv.e2=glv.e^2;
glv.Rp=(1-glv.f)*glv.Re;%�̰���
glv.ep=sqrt(glv.Re^2-glv.Rp^2)/glv.Rp;%�ڶ�ƫ����,  �˴��иĶ����ӺŸ�Ϊ����
glv.ep2=glv.ep^2;
glv.wie=7.2921151467e-5;%������ת������
glv.g0=9.7803267714;%�������ٶ�
glv.mg=1.0e-3*glv.g0;%���������ٶ�
glv.ug=1.0e-6*glv.g0;%΢�������ٶ�
glv.ppm=1.0e-6;%�����֮һ
glv.deg=pi/180;%�Ƕ�
glv.min=glv.deg/60;%�Ƿ�
glv.sec=glv.min/60;%����
glv.hur=3600;%Сʱ
glv.dph=glv.deg/glv.hur;%��/Сʱ
glv.cs=[    %Բ׶�˶��ͻ���ЧӦ����ϵ��
    2./3,      0,         0,         0
    9./20,     27./20,    0,         0
    54./105,   92./105,   214./105,  0
    250./504,  525./504,  650./504,  1375./504];