function [ts,T_diff_r,T_diff_R,ts90_r,ts90_R,A,x,y]=get_transient_time_analytical(alpha,phi_o,r_ratio,KsMr,KfMr,KlMr,MmMr,L2kappa,EtaMr)

%This function is used to calculate the postinjection evolution time scale
%ts. Note that when running this function, plot(x,y) each time and make sure that it
%is indeed the smallest x-value that gets picked up. 

%This code is used to calculate the results shown in Figure 5 and Figure8.

%L2kappa: either r_o^2/kappa (when r_o is fixed), or R_o^2/kappa (when R_o
%is fixed) 
%EtaMr: eta_f/\mu_r, ratio between melt viscosity (Pa.s) and rock rigidity
%(Pa). This term has time-scale (second).  


KsMm=KsMr/MmMr;
MrMm=1/MmMr;
KfKs=KfMr/KsMr;

KuMm=(1-alpha)*KsMm+alpha^2*KfKs*KsMm/(phi_o+(alpha-phi_o)*KfKs);
KmMm=(1-alpha)*KsMm;
KmMr=KmMm*MmMr;
KuMr=KuMm*MmMr;

KfMm=KfKs*KsMm;
a1=(KmMm+4/3)*(KuMm-KmMm)/(alpha^2*MrMm*(KuMm+4/3));
a2=-(KuMm-KmMm)/(alpha*(KuMm+4/3));
b1=-(r_ratio^(2/3))*(KuMm-KmMm)/(alpha*(KuMm+4/3));
b2=-(r_ratio^(-1/3))*(r_ratio-1)*MrMm/(3*(KuMm+4/3));
c1=4*r_ratio*MmMr*(KuMm-KmMm)/(alpha*(KuMm+4/3));
c2=1+(r_ratio-1)*(4/3)/(KuMm+4/3);
d1=-3*r_ratio*KlMr-b1*3*KlMr*(r_ratio^(1/3));
d2=-b2*3*KlMr*(r_ratio^(1/3));
e0=4*(MmMr-1);
f1=a1;
 f2=a2*e0-d2*e0+3*r_ratio*KlMr;
 f3=-d1;
% f4=delta*KlMr;   %%%%%linear: Delta-dependent
 g1=c1+d1;
 g2=c2*e0+d2*e0-4*MmMr*r_ratio-3*r_ratio*KlMr;
h0=f3-f2*g1/g2;
%h1=f4+f2*f4/g2;    %%%%linear: Delta-dependent

r_o=r_ratio^(-1/3);
x=0:0.1:100;
y=sin(sqrt(x)*(1-r_o)).*((f1/r_o-h0*r_o)*x-h0)-sqrt(x).*(f1*(x)/r_o+h0*r_o-h0).*cos(sqrt(x)*(1-r_o));

%fun=@(x) Laplace_root_2(x,r_o,h0,f1);   %when the shell is too thin
fun=@(x) Laplace_root_2(1/x,r_o,h0,f1);
%ts=fzero(fun,0.5); %this will give the largest value of 1/x, for root x
ts=fzero(fun,1/10);
%ts=1/x;
%!! Only if time scale is smaller than 0.5\tau_diff
ts90=-log(0.1)*ts;
T_diff_r=L2kappa*EtaMr*r_ratio^(2/3)*alpha^2*(KuMr+4*MmMr)/((KmMr+4*MmMr/3)*(KuMr-KmMr));
T_diff_R=L2kappa*EtaMr*alpha^2*(KuMr+4*MmMr)/((KmMr+4*MmMr/3)*(KuMr-KmMr));
ts90_r=ts90*T_diff_r;
ts90_R=ts90*T_diff_R;

%mush-dependent factor in the diffusion time-scale
A=alpha^2*(KuMr+4*MmMr/3)/((KmMr+4*MmMr/3)*(KuMr-KmMr));


%a=phi_o*(KmMr+4*MmMr/3)*(KuMr-KmMr)/(alpha^2*KfMr*(KuMr+4*MmMr/3));
% c_{diff}=a*c_{Darcy}
