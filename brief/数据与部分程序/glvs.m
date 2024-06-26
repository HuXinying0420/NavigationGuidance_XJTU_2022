%全局变量
global glv
glv.Re=6378160;%地球半径(长半轴)
glv.f=1/298.3;%地球扁率
glv.e=sqrt(2*glv.f-glv.f^2);%偏心率
glv.e2=glv.e^2;
glv.Rp=(1-glv.f)*glv.Re;%短半轴
glv.ep=sqrt(glv.Re^2-glv.Rp^2)/glv.Rp;%第二偏心率,  此处有改动，加号改为减号
glv.ep2=glv.ep^2;
glv.wie=7.2921151467e-5;%地球自转角速率
glv.g0=9.7803267714;%重力加速度
glv.mg=1.0e-3*glv.g0;%毫重力加速度
glv.ug=1.0e-6*glv.g0;%微重力加速度
glv.ppm=1.0e-6;%百万分之一
glv.deg=pi/180;%角度
glv.min=glv.deg/60;%角分
glv.sec=glv.min/60;%角秒
glv.hur=3600;%小时
glv.dph=glv.deg/glv.hur;%度/小时
glv.cs=[    %圆锥运动和划桨效应补偿系数
    2./3,      0,         0,         0
    9./20,     27./20,    0,         0
    54./105,   92./105,   214./105,  0
    250./504,  525./504,  650./504,  1375./504];