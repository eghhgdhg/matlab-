 alphar=(0:10:710);
 p1g=[0.0857,0.08941,0.08336,0.08385,0.08229,0.0818,0.08414,0.08453, 0.08346,0.0858,0.08629,0.08751,0.08921,0.09195,0.09239,0.09292, 0.09478,0.09582,0.098,0.0986,0.09982,0.10336,0.10736,0.11673,  0.13087,0.15526,0.17838,0.21916,0.27409,0.35965,0.49906,0.72989,  1.1332,1.87837,3.14011,4.69971,6.70251,6.73661,5.11565,3.51069,   2.41675,1.67978,1.2474,0.95628,0.76555,0.64058,0.55072,0.48428, 0.44136,0.40536,0.38009,0.36594,0.34233,0.31502,0.20931,0.15575, 0.1379,0.11439,0.09486,0.08219,0.09156,0.08953,0.08882,0.09741, 0.10317,0.09468,0.08941,0.08424,0.08541,0.09302,0.09673,0.0936];
 m3=1.38;n=2000;m4=1.81;d=0.095;mk1=0.563;mk2=1.15;rho=0.0415;r=0.0575;l=0.175;m1=0.543;P0=0.101325; 
 w=2*pi*n/60;
	  lambda=r/l;    % 活塞组质量m3 连杆组质量m4  曲柄销质量mk1 曲柄臂质量mk2  连杆小头质量M1[i] 连杆大头质量m2  旋转质量mr  往复运动质量mj
	  m2=m4-m1;
	  mk=mk1+2*rho/r*mk2;     %曲柄臂质心回转半径rho  曲柄半径r  连杆长度l 
	  mr=m2+mk;
	  mj=m1+m3;
	  A=pi/4*d*d;

     alpha=alphar*pi/180; 
	 beta=asin(lambda*sin(alpha)) ;     
	 bet=180*beta/pi;
	 B= 1-lambda^2*(sin(alpha)).^2;
      wl=lambda*w*cos(alpha)./sqrt(B); 
      al=w*w*lambda*(lambda^2-1)*sin(alpha)./B.^(3/2);      %角速度w(rad/s)  曲轴转角alpha   连杆运动角速度wl   
      
    x=r*(1/lambda+1-cos(alpha)-1/lambda*cos(beta));    
 	v=r*w*(sin(alpha+beta)./cos(beta));
 	a=w*w*r*(cos(alpha+beta))./cos(beta)+lambda*(cos(alpha)).^2./(cos(beta)).^3;
     pj=-a*mj/A*10^(-6);   
    pg=p1g-P0; 
    p=pg+pj;                 %气缸上部绝对压力p1g    大气压力p0  相对压力pg   活塞往复惯性压力pj   活塞上总压力p  侧压力pn  连杆力pl   切向力t 径向力k 
    pn=p.*tan(beta);
    pl=p./cos(beta);
    t=pl.*sin(alpha+beta);
    k=pl.*cos(alpha+beta);
    M1=t*r*10^6*A;                                         
	Mz1=0;
	Mz2=M1;
	Mq1=Mz1+0.5*M1;	
    
 	krl=m2*r*w*w/A*10^(-6);				%连杆旋转质量m2的离心力krl   连杆作用于曲柄销上的负荷pq 	alpha q1  pq对应的角度 alpha q   连杆轴承负荷pp  	连杆轴承负荷的对应的角度alphap 							   
	pqx=krl-k;
	pqy=t;													   																									   
	pq=sqrt(pqx.*pqx+pqy.*pqy);
	alphaq1=atan(abs(pqy./pqx));
	alphaq=(1-sign(pqx))*pi/2+sign(pqx)*sign(pqy).*alphaq1;		
	pp=pq;
	alphap=alphaq+pi+alpha+beta;									           
	
	krk=mk*r*w*w/A*10^(-6);     %mk旋转产生的离心力产生的压强krk    主轴颈上的负荷pz    主轴颈上的负荷对应的角度alpha z  pzx  pzy   主轴承上的负荷pc   主轴承上的负荷对应的角度alpha c 
	pzx=-0.5*(krl+krk-k);
	pzy=-0.5*t;
	pz=sqrt(pzx.*pzx+pzy.*pzy);
	alphaz1=atan(abs(pzy./pzx));
	alphaz=(1-sign(pzx))*pi/2+sign(pzx)*sign(pzy)*alphaz1;
	pc=pz;
	alphac=alphaz+pi+alpha;
 %       polarplot(alphaz,pz,'r-','linewidth',2);
%        title('主轴颈负荷图');
             polarplot(alphac,pc,'b-','linewidth',2);
 title('主轴承负荷图');
  %   polarplot(alphaq,pq,'r-','linewidth',2);
 %title('连杆负荷图');
 %thetalim([90 150]);
% thetaticks(90:5:150);
  %      polarplot(alphap,pp,'r-','linewidth',2);
 %title('连杆轴承负荷图');
    
    function   z= sign(x)
 z=x/abs(x);
    end
 