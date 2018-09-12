function [S1,u1,u0,Pl,u1_vis_fin,u0_vis_fin,Pl_vis_fin]=get_viscoelastic(t,alpha,phi_o,r_ratio,KsMr,KfMr,KlMr,MmMr,delta)
% This code is used to solve for the evoltion of a viscoelastic shell
% (undrained poroviscoelastic) after a sudden injection. 
%In the input:

% r_ratio=R_o^3/r_o^3; 
% KsMr is the ratio between solid bulk modulus K_s and regidity \Mu_r
% KfMr is the ratio between pore fluid bulk modulus K_f and regidity \Mu_r
% KlMr is the ratio between core fluid bulk modulus K_l and regidity \Mu_r
% MmMr is the ratio between mush regidity \mu_m and regidity \Mu_
% delta=dM/Mo is the ratio between injected magma and original core magma
% All results are normalized by their respective scales
% Note that time is normalized by the viscoelastic relaxation time
% \eta_m/\mu_m. This is why we do not input a value eta_m, as it only
% determined the time scale. 

%Relevant output:
%Pl: core pressure as function of time 
%u1: displacement at r=R_o, which is proportional to tensile stress
% S1: the post-injection time scale (normalized)

%This code was used to compute results shown in Figure 8 


%According to the alaytical solution, the solution u(1,t)=A0+A1exp(-S1*t), t normalized by \eta_m/\mu_m


KsMm=KsMr/MmMr;
KlMm=KlMr/MmMr;
MrMm=1/MmMr;
KfKs=KfMr/KsMr;

KuMm=(1-alpha)*KsMm+alpha^2*KfKs*KsMm/(phi_o+(alpha-phi_o)*KfKs);

KmMm=(1-alpha)*KsMm;
KfMm=KfKs*KsMm;
KuMr=KuMm*MmMr;

A=4*(KuMr+KlMr*(r_ratio-1))+3*KlMr*KuMr*r_ratio;
B=4*MmMr*(KuMr-KlMr-KuMr*r_ratio)-16*MmMr*r_ratio/3-4*(KuMr+KlMr*(r_ratio-1))-3*KlMr*KuMr*r_ratio;

S1=A/B;  % this corresponds to the viscoelastic time scale \eta_m/\mu_m
%S2=-KuMr/(KuMr-4*MmMr/3);

I1=4*(KuMr+KlMr*(r_ratio-1))+3*KlMr*KuMr*r_ratio;  %I1 is normalized by mu_r
I0=4*(KuMr-KlMr-KuMr*r_ratio)-16*r_ratio/3;

A0_u=KlMr*delta*KuMr/I1; 
A1_u=-KlMr*delta*((KuMr+4*MmMr/3)/(MmMr*I0-I1)+KuMr/I1);
A0_Phi=1/I1;
A1_Phi=(MmMr-1)/(I0*MmMr-I1)-1/I1;
% A1_Phi=((KuMr+4*MmMr/3)/(MmMr*I0-I1)+KuMr/I1)...
%     *(I1*MmMr*(MmMr-1)-MmMr*(I0*MmMr-I1))/(I1*MmMr*(KuMr-4*MmMr/3)+KuMr*MmMr*(I0*MmMr-I1));
% A2_Phi=(8/3)*MmMr*KuMm...
%     *(1/KuMr+(MmMr-1)/(KuMr-4*MmMr/3))/(I1*(KuMr-4*MmMr/3)+KuMr*(I0*MmMr-I1));
% 
u1=A0_u+A1_u*exp(S1.*t);
u1_vis_fin=A0_u;

Phi=A0_Phi+A1_Phi*exp(S1.*t);
Phi_fin=A0_Phi;
%Phi=A0_Phi+A1_Phi*exp(S1.*t)+A2_Phi*exp(S2.*t);

u0=(4*KlMr*delta/3)*(r_ratio^(2/3)-r_ratio^(-1/3))*Phi+r_ratio^(2/3)*u1;
u0_vis_fin=(4*KlMr*delta/3)*(r_ratio^(2/3)-r_ratio^(-1/3))*Phi_fin+r_ratio^(2/3)*u1_vis_fin;

Pl=-4*KlMr*KlMr*delta*(r_ratio-1)*Phi-3*KlMr*r_ratio*u1+KlMr*delta;
Pl_vis_fin=-4*KlMr*KlMr*delta*(r_ratio-1)*Phi_fin-3*KlMr*r_ratio*u1_vis_fin+KlMr*delta;


