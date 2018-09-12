 function [M_s,u1_fin,u1_ini,u0_ini,u0_fin,tensial_fin,tensial_ini,Pl_fin,Pl_ini,Pf_fin,Pf_ini,P_ini,P_fin,phi_ini,phi_fin]=IniFin_values2(alpha,phi_o,r_ratio,KsMr,KfMr,KlMr,MmMr,delta);
% This code calcualte the instantaneous response of the system to a sudden
% injection (initial), as well as the final, steady state response after
% the system has evolved enough time (final)

%In the input,
% r_ratio=R_o^3/r_o^3; 
% KsMr is the ratio between solid bulk modulus K_s and regidity \Mu_r
% KfMr is the ratio between pore fluid bulk modulus K_f and regidity \Mu_r
% KlMr is the ratio between core fluid bulk modulus K_l and regidity \Mu_r
% MmMr is the ratio between mush regidity \mu_m and regidity \Mu_
% delta=dM/Mo is the ratio between injected magma and original core magma
% All results are normalized by their respective scales

%In the output:
% *_ini indicates the instantaneous response at t=0; 
% *_fin indicates the final response at t=infinity;
% The code is based on the analytical solution using Laplace transform.
% The code was used to compute results shown in Figure 6,7 and FigureA3,A4


 
 
%%% Similar to IniFin_values but all stresses normalized by \mu_r instead of
%%% \mu_m
 %[a1,a2,b1,b2,c1,c2,d1,d2,e0,f1,f2,f3,h0,h1]..
%See Notes
% alpha;
%r_ratio=R^3/r^3;
%KsMm: K_s/mu_m
%KfKs: K_f/K_s
%KlMm: K_l/mu_m
%delta:dM/Mo
%MrMm: mu_r/mu_m

%Pf=Pf/mu_m;
%Pl=Pl/mu_m;
%drhol=d \rho/\rho_o in the core fluid
%drhof=d \rho/\rho_o in the pore fluid
%dVcore=dVcore/Vcore_o;
%dVtot=dVtot/Vtot_o;
%tensial=tensial stress on rock (u=Ro)/mu_m





KsMm=KsMr/MmMr;
KlMm=KlMr/MmMr;
MrMm=1/MmMr;
KfKs=KfMr/KsMr;

KuMm=(1-alpha)*KsMm+alpha^2*KfKs*KsMm/(phi_o+(alpha-phi_o)*KfKs);

KmMm=(1-alpha)*KsMm;
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
f4=delta*KlMr;   %%%%%linear: Delta-dependent
g1=c1+d1;
g2=c2*e0+d2*e0-4*MmMr*r_ratio-3*r_ratio*KlMr;
h0=f3-f2*g1/g2;
h1=f4+f2*f4/g2;    %%%%linear: Delta-dependent

m0=h1/(f1+h0*(1-1/r_ratio)/3);    %%%linear: Delta-dependent

u1_ini=-f4/g2; %%%linear: Delta-dependent
u1_fin=-f4/g2-g1*(h1/h0-f1*m0/h0)/g2;    %%%linear: Delta-dependent

zeta_ini=-4*f4*(MmMr-1)/g2;     %%%linear: Delta-dependent
zeta_fin=4*(MmMr-1)*u1_fin;   %%%linear: Delta-dependent

u0_ini=b2*zeta_ini+r_ratio^(2/3)*u1_ini;
u0_fin=b1*(1-1/r_ratio)*m0/3+b2*zeta_fin+r_ratio^(2/3)*u1_fin;

Pl_ini=-3*r_ratio*KlMr*u1_ini+d2*zeta_ini+KlMr*delta;   %%%linear: Delta-dependent
Pl_fin=d1*(1-1/r_ratio)*m0/3+d2*zeta_fin-3*KlMr*r_ratio*u1_fin+KlMr*delta;   %%%linear: Delta-dependent

Pf_ini=a2*zeta_ini;    %%%linear: Delta-dependent
Pf_fin=a1*m0+a2*zeta_fin;    %%%linear: Delta-dependent

%m_ini_edge=f4/f1-f2*u1_ini/f1;
%Pf_ini_edge=a1*m_ini_edge+a2*zeta_ini;
%dP_ini_alt=(d2-a2)*zeta_ini-3*KlMm*r_ratio*u1_ini+KlMm*delta;

P_ini=-zeta_ini*KuMm/(KuMm+4/3);    %%%linear: Delta-dependent
P_fin=(4/3)*MmMr*(KuMm-KmMm)*m0/(alpha*(KuMm+4/3))-zeta_fin*KuMm/(KuMm+4/3);   %%%linear: Delta-dependent

tensial_ini=2*u1_ini;    %%%linear: Delta-dependent
tensial_fin=2*u1_fin;    %%%linear: Delta-dependent

% phi_ini=-phi_o*Pf_ini/KfMm;    %%%linear: Delta-dependent
% phi_fin=m0-phi_o*Pf_fin/KfMm;   %%%linear: Delta-dependent
phi_ini=(1-(1-phi_o)/(1-alpha))*(P_ini-Pf_ini)/KsMr;
phi_fin=(1-(1-phi_o)/(1-alpha))*(P_fin-Pf_fin)/KsMr;



M_s=3*r_ratio*h1*(1-1/r_ratio)./(3*f1+h0*(1-1/r_ratio));    %%%linear: Delta-dependent

if (r_ratio==1)   %Pure elastic 
    M_s=0;
    Pf_fin=NaN;
    Pf_ini=NaN;
    P_ini=NaN;
    P_fin=NaN;
    phi_ini=NaN;
    phi_fin=NaN;
end
