clear;clc;
tic
%[Lp Wp Dr L W Lw Lf Lr Rmin Omax dOmax]
%(7.5,2.5,6,4.718,1.84,2.74,0.968,1.01,6,0.428,0.43)
[Lp,Wp,Dr,L,W,Lw,Lf,Lr,Rmin,Omax,dOmax]=...
    deal(7.5,2.5,6,4.718,1.84,2.74,0.968,1.01,6,0.428,0.43);
[x1,y1,Q1,O1,x3,y3,Q3,O3]=...
    deal(-6,2.75,0,0,-(Lp/2-Lr-0.2),0,0,Omax);

[Q2,O2]=deal(0,0);

syms x xx4 yy4 QQ4
f7=xx4+Lw*sin(QQ4)-W*cos(QQ4)/2-Wp/2;
f8=x3+Rmin*sin(QQ4)-xx4;
f9=y3+Rmin*(1-cos(QQ4))-yy4;
[x4_1,y4_1,Q4_1]=solve(f7,f8,f9,xx4,yy4,QQ4);
x4=double(x4_1(Q4_1<pi/2));y4=double(y4_1(Q4_1<pi/2));
Q4=double(Q4_1(Q4_1<pi/2));O4=O3;
Y=Rmin-(Rmin^2-(x-x3)^2)^(1/2);dY=diff(Y);
syms  a0 a1 a2 a3 a4 a5
num=[a0 a1 a2 a3 a4 a5];
y=poly2sym(num);dy=diff(y);ddy=diff(dy);
f1=y1-subs(y,x1);
f3=tan(Q1)-subs(dy,x1);
f5=Lw*subs(ddy,x1)-(1+(subs(dy,x1))^2)^(3/2)*tan(O1);
f12=y4-subs(y,x4);
f14=tan(Q4)-subs(dy,x4);
f16=Lw*subs(ddy,x4)-(1+(subs(dy,x4))^2)^(3/2)*tan(O4);
%以下需使用循环
xlow=Lp/2+4.5;xup=Lp*2;ylow=Wp/2+W/2+1;yup=Wp/2+Dr-W/2;
xdim=5;ydim=5;
xm=linspace(xlow,xup,xdim);ym=linspace(ylow,yup,ydim);
xmm=linspace(1,1,length(ym))'*xm;
ymm=ym'*linspace(1,1,length(xm));
sp(:,1)=(xmm(:));
sp(:,2)=(ymm(:));

i=1;
while((i>0)&&(i<length(sp)))
    fprintf('%d ',i);
    [x2,y2]=deal(sp(i,1),sp(i,2));
    f2=y2-subs(y,x2);
    f4=tan(Q2)-subs(dy,x2);
    f6=Lw*subs(ddy,x2)-(1+(subs(dy,x2))^2)^(3/2)*tan(O2);
    x_1=solve(f1,f2,f3,f4,f5,f6,a0,a1,a2,a3,a4,a5);
    f11=y2-subs(y,x2);
    f13=tan(Q2)-subs(dy,x2);
    f15=Lw*subs(ddy,x2)-(1+(subs(dy,x2))^2)^(3/2)*tan(O2);
    x_2=solve(f11,f12,f13,f14,f15,f16,a0,a1,a2,a3,a4,a5);
    %离散化
    inter=20;
    x_p1=linspace (x1,x2,inter); x_p3=linspace(x4,x2,inter);
    if(any(abs(subs(subs(ddy/(1+dy^2)^(3/2),x_1),x_p1))>(1/Rmin)))
        i=i+1;
        continue  
    else
        if(any(abs(subs(subs(ddy/(1+dy^2)^(3/2),x_2),x_p3))>(1/Rmin)))
            i=i+1;
            continue
        else
            y_p1=subs(subs(y,x_1),x_p1);Q_p1=atan(double(subs(subs(dy,x_1),x_p1)));
            y_A_p1=y_p1+(Lf+Lw)*sin(Q_p1)+W*cos(Q_p1)/2;y_C_p1=y_p1-Lr*sin(Q_p1)-W*cos(Q_p1)/2;
            y_B_p1=y_p1+(Lf+Lw)*sin(Q_p1)-W*cos(Q_p1)/2;y_D_p1=y_p1-Lr*sin(Q_p1)+W*cos(Q_p1)/2;
            if((any(y_A_p1>Wp/2+Dr)||any(y_B_p1<Wp/2)||any(y_C_p1<Wp/2)||any(y_D_p1>Wp/2+Dr)))
                i=i+1;
                continue
            else
                y_p3=subs(subs(y,x_2),x_p3);Q_p3=atan(double(subs(subs(dy,x_2),x_p3)));
                y_A_p3=y_p3+(Lf+Lw)*sin(Q_p3)+W*cos(Q_p3)/2;y_C_p3=y_p3-Lr*sin(Q_p3)-W*cos(Q_p3)/2;
                y_B_p3=y_p3+(Lf+Lw)*sin(Q_p3)-W*cos(Q_p3)/2;y_D_p3=y_p3-Lr*sin(Q_p3)+W*cos(Q_p3)/2;
                if((any(y_A_p3>Wp/2+Dr)||any(y_D_p3>Wp/2+Dr)))
                    i=i+1;
                    continue
                else
                    x_p2=linspace(x3,x4,inter);y_p2=subs(Y,x_p2);Q_p2=atan(double(subs(dY,x_p2)));
                    y_A_p2=y_p2+(Lf+Lw)*sin(Q_p2)+W*cos(Q_p2)/2;y_C_p2=y_p2-Lr*sin(Q_p2)-W*cos(Q_p2)/2;
                    y_B_p2=y_p2+(Lf+Lw)*sin(Q_p2)-W*cos(Q_p2)/2;y_D_p2=y_p2-Lr*sin(Q_p2)+W*cos(Q_p2)/2;
                    x_t=[x_p1 fliplr(x_p3) fliplr(x_p2)];
                    y_t=[y_p1 fliplr(y_p3) fliplr(y_p2)];
                    Q_t=[Q_p1 fliplr(Q_p3) fliplr(Q_p2)];
                    x_A=x_t+(Lf+Lw)*cos(Q_t)-W*sin(Q_t)/2;y_A=[y_A_p1 fliplr(y_A_p3) fliplr(y_A_p2)];
                    x_B=x_t+(Lf+Lw)*cos(Q_t)+W*sin(Q_t)/2;y_B=[y_B_p1 fliplr(y_B_p3) fliplr(y_B_p2)];
                    x_C=x_t-Lr*cos(Q_t)+W*sin(Q_t)/2;y_C=[y_C_p1 fliplr(y_C_p3) fliplr(y_C_p2)];
                    x_D=x_t-Lr*cos(Q_t)-W*sin(Q_t)/2;y_D=[y_D_p1 fliplr(y_D_p3) fliplr(y_D_p2)];
                    toc
                    figure(1);hold on;
                    %             plot(sp(:,1),sp(:,2),'o');
                    %             plot(sp(:,1),sp(:,2));
                    %             plot(x_t,y_t);
                    plot([x1-3 xup],[Wp/2 Wp/2],'b');
                    plot([x1-3 xup],[Wp/2+Dr/2 Wp/2+Dr/2],'b','linestyle','--');
                    plot([x1-3 xup],[Wp/2+Dr Wp/2+Dr],'b');
                    plot([-Lp/2 -Lp/2],[-Wp/2 Wp/2],'b');
                    plot([Lp/2 Lp/2],[-Wp/2 Wp/2],'b');
                    plot([-Lp/2 Lp/2],[-Wp/2 -Wp/2],'b');
                    for ii=1:inter*3
                        plot([x_A(ii) x_B(ii)],[y_A(ii) y_B(ii)],'k');
                        plot([x_B(ii) x_C(ii)],[y_B(ii) y_C(ii)],'k');
                        plot([x_C(ii) x_D(ii)],[y_C(ii) y_D(ii)],'k');
                        plot([x_D(ii) x_A(ii)],[y_D(ii) y_A(ii)],'k');
                    end
                    axis equal;
                    break
                end
            end
        end
    end
end














