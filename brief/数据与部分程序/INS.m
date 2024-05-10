clc;clear;close all;
%�����ߵ����·��棬��Ԫ������һ�׽����㷨��������Բ׶�����ͻ�������


ts=0.1;%����ʱ��

Re=6378160;%���򳤰���
wie=7.2921151467e-5;%������ת������
f=1/298.3;%�������
g0=9.7803;%�������ٶȻ���ֵ
deg=pi/180;%�Ƕ�
min=deg/60;%�Ƿ�
sec=min/60;%����
hur=3600;%Сʱ
dph=deg/hur;%��/Сʱ

%��ȡ����
wbibS=dlmread('dataWbibN.txt');
fbS=dlmread('dataFbibN.txt');
posS=dlmread('dataPos.txt');
vtetS=dlmread('dataVn.txt');
p_gps=dlmread('dataGPSposN.txt');

%ͳ�ƾ����ʼ��
[mm,nn]=size(posS);
posSta=zeros(mm,nn);
vtSta=posSta;
attSta=posSta;

posC(:,1)=posS(:,1);  %λ��������ʼֵ
vtC(:,1)=vtetS(:,1);  %�ٶ�������ʼֵ
attC(:,1)=[ 0;
                0;
          0.3491];  %��̬��������ʼֵ


Qk=1e-6*diag([0.01,0.01,0.01,0.01,0.01,0.01,0.9780,0.9780,0.9780]).^2;%ϵͳ�����������
Rk=diag([1e-5,1e-5,10.3986]).^2;   %������������

Tg = 3600*ones(3,1);             %������Markov�������ʱ��
Ta = 1800*ones(3,1);             %���ٶȼ�Markov�������ʱ��

GPS_Sample_Rate=10;  %GPS������̫��Ч��Ҳ����


StaNum=10;%�ظ����д�����������ȡͳ��ƽ��ֵ

for OutLoop=1:StaNum
    
    Pk = diag([min,min,min, 0.5,0.5,0.5, 30./Re,30./Re,30, 0.1*dph,0.1*dph,0.1*dph, 0.1*dph,0.1*dph,0.1*dph, 1.e-3,1.e-3,1.e-3].^2); %��ʼ����Э�������           
    Xk = zeros(18,1);  %��ʼ״̬
    %%
    N=size(posS,2);
%     j=1; 
    for k=1:N-1
        si=sin(attC(1,k));ci=cos(attC(1,k));
        sj=sin(attC(2,k));cj=cos(attC(2,k));
        sk=sin(attC(3,k));ck=cos(attC(3,k));
        %kʱ����̬����
        M=[cj*ck+si*sj*sk,    ci*sk,  sj*ck-si*cj*sk;
            -cj*sk+si*sj*ck,  ci*ck,  -sj*sk-si*cj*ck;
            -ci*sj,           si,     ci*cj]; %��Cnb����
        q_1=[                     sqrt(abs(1.0+M(1,1)+M(2,2)+M(3,3)))/2.0;
            sign(M(3,2)-M(2,3))*sqrt(abs(1.0+M(1,1)-M(2,2)-M(3,3)))/2.0;
            sign(M(1,3)-M(3,1))*sqrt(abs(1.0-M(1,1)+M(2,2)-M(3,3)))/2.0;
            sign(M(2,1)-M(1,2))*sqrt(abs(1.0-M(1,1)-M(2,2)+M(3,3)))/2.0];
        fn(:,k)=M*fbS(:,k);%����������任
        
       
        
        %�����ߵ�����
        wnie=wie*[0;cos(posC(1,k));sin(posC(1,k))];%����ϵ��Թ���ϵ��ת�����ٶ��ڵ���ϵ������ϵ����ͶӰ
        %���������ʰ뾶
        Rm=Re*(1-2*f+3*f*sin(posC(1,k))^2)+posC(3,k);
        Rn=Re*(1+f*sin(posC(1,k))^2)+posC(3,k);

        wnen=[-vtC(2,k)/(Rm+posC(3,k));vtC(1,k)/(Rn+posC(3,k));vtC(1,k)*tan(posC(1,k))/(Rn+posC(3,k))];%����ϵ�����Ե���ϵ��ת�����ٶ��ڵ���ϵ��ͶӰ
        g=g0+0.051799*sin(posC(1,k))^2-0.94114e-6*posC(3,k);%�������ٶ�
        gn=[0;0;-g];%��������ϵ���������ٶ�
       
       
        wbnbC(:,k)=wbibS(:,k)-M\(wnie+wnen); %��̬��ת�������ʼ���
        q=1.0/2*qmul(q_1,[0;wbnbC(:,k)])*ts+q_1;  %��̬��Ԫ������
        q=q/sqrt(q'*q);%��Ԫ����һ��
        
        %��̬�Ǹ���
        q11=q(1)*q(1);q12=q(1)*q(2);q13=q(1)*q(3);q14=q(1)*q(4);
        q22=q(2)*q(2);q23=q(2)*q(3);q24=q(2)*q(4);
        q33=q(3)*q(3);q34=q(3)*q(4);
        q44=q(4)*q(4);
        T=[q11+q22-q33-q44,  2*(q23-q14),      2*(q24+q13);
            2*(q23+q14),     q11-q22+q33-q44,  2*(q34-q12);
            2*(q24-q13),     2*(q34+q12),      q11-q22-q33+q44];
    
        attC(:,k+1)=[asin(T(3,2));atan(-T(3,1)/T(3,3));atan(T(1,2)/T(2,2))];
        %�����gamma����
        if T(3,3)<0
           if attC(2,k+1)<0
               attC(2,k+1)=attC(2,k+1)+pi;
           else
               attC(2,k+1)=attC(2,k+1)-pi;
           end
        end
        %�����psi����
        if T(2,2)<0
           if T(1,2)>0
               attC(3,k+1)=attC(3,k+1)+pi;
           else
               attC(3,k+1)=attC(3,k+1)-pi;
           end
        end
        if abs(T(2,2))<1.0e-20
           if T(1,2)>0
              attC(3,k+1)=pi/2.0;
           else
              attC(3,k+1)=-pi/2.0;
           end
        end
        
        %�ٶȸ���
        vtC(:,k+1)=vtC(:,k)+ts*(fn(:,k)+gn-cross((2*wnie+wnen),vtC(:,k)));
        %λ�ø���
        posC(1,k+1)=posC(1,k)+ts*vtC(2,k)/(Rm+posC(3,k));                 %γ��
        posC(2,k+1)=posC(2,k)+ts*vtC(1,k)/((Rn+posC(3,k))*cos(posC(1,k)));%����
        posC(3,k+1)=posC(3,k)+ts*vtC(3,k);        %�߶�
        
        %��ӳ���
        [ F,G,H ] = GetConSis( vtC(:,k), posC(:,k), q , fn(:,k), Tg, Ta );%����GetConSis�õ�����ϵͳ״̬ת�ƾ���F������ϵͳ��������G������ϵͳ��������K���� 
        %ins����ٶ�����\λ������\��̬��Ԫ��\�������룬����ϵ�±�������\���ٶȼ����Ư�����ʱ��\���������Ư�����ʱ��
        
        if mod(k+1,10) == 0 %�ж�k�Ƿ�Ϊ10�ı������������Kalman�˲�����������
            Zk = - posC(:,k+1) + p_gps(:,(k+1)/10);%Zk  ����λ���������
            [ E_att, E_pos, E_vn, Xk, Pk ] = kalman_GPS_INS_correct( Xk, Qk, Pk, F, G, H, ts , Zk , Rk);%Kalman�˲�
            %k-1ʱ��״̬����������ֵ��\ϵͳ�����������\���ƾ���������  ��������\����λ���������\���������������
            
            %����λ�á��ٶ�����
            posC(:,k+1) = posC(:,k+1) + E_pos;
            vtC(:,k+1) = vtC(:,k+1) + E_vn;
        else
            [ E_att, E_pos, E_vn, Xk, Pk ] = kalman_GPS_INS_correct( Xk, Qk, Pk, F, G, H, ts );%GPS�������Ե���
        end

    end
    

        
    %ͳ�ƾ������
    attSta=attSta+attC;
    vtSta=vtSta+vtC;
    posSta=posSta+posC;
    
end
%��ͳ�ƾ���ȡƽ��
attC=1./StaNum*attSta;
posC=1./StaNum*posSta;
vtC=1./StaNum*vtSta;

%����ֵ�����ֵ����� ��ͼ
%��������Ϊ3x(N+1)����������
%N�����ݣ���������������Ϊ0.1s,����ͼ��������
for i=1:N
    time(i)=i*ts;
    Rm = Re*(1-2*f+3*f*(sin(attC(1,i)))^2);
    Rn = Re*(1+f*(sin(attC(1,i)))^2);
    Latitude_error(i)=(posC(1,i)-posS(1,i))*Rm;
    Longitude_error(i)=(posC(2,i)-posS(2,i))*Rn*cos(attC(1,i));
end

posCp=posC(:,1:N);
figure;
subplot(131);plot(time,Latitude_error);title('γ�����');xlabel('Time /s');ylabel('\it\deltaL /m');grid on;
subplot(132);plot(time,Longitude_error);title('�������');xlabel('Time /s');ylabel('\it\delta\phi /m');grid on;
subplot(133);plot(time,posCp(3,:)-posS(3,:));title('�߶����');xlabel('Time /s');ylabel('\it\deltah /m');grid on;

vtCp=vtC(:,1:N);
figure;
subplot(131);plot(time,vtCp(1,:)-vtetS(1,:));title('���ٶ����');xlabel('Time /s');ylabel('\it\deltavelocity east /(m/s)');grid on;
subplot(132);plot(time,vtCp(2,:)-vtetS(2,:));title('���ٶ����');xlabel('Time /s');ylabel('\it\deltavelocity north /(m/s)');grid on;
subplot(133);plot(time,vtCp(3,:)-vtetS(3,:));title('���ٶ����');xlabel('Time /s');ylabel('\it\deltavelocity up /(m/s)');grid on;

%��ά���й켣ͼ
figure;
plot3(posS(2,:)/pi*180,posS(1,:)/pi*180,posS(3,:),'k');
hold on;
plot3(posCp(2,:)/pi*180,posCp(1,:)/pi*180,posCp(3,:),'r');grid on;
ylabel('γ��L /arcdeg');xlabel('����\phi /arcdeg');zlabel('�߶�h /m');title('����-���������й켣������-INS������й켣');
