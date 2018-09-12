
function [S1,u1,u0]=get_viscoelastic_gradual(t,alpha,phi_o,r_ratio,KsMr,KfMr,KlMr,MmMr,delta,tau)
%Here the time is normalized by the  relaxation time 
%the solution u(1,t)=A0+A1exp(-S1*t), t normalized by the relaxation time \eta_m/\mu_m


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

I1=4*(KuMr+KlMr*(r_ratio-1))+3*KlMr*KuMr*r_ratio;  %I1 is normalized by mu_r
I0=4*(KuMr-KlMr-KuMr*r_ratio)-16*r_ratio/3;

S1=A/B;  % this corresponds to the viscoelastic time scale \eta_m/\mu_m
X1=-S1;

D=-exp(-X1*(t-tau).*heaviside(t-tau))/X1^2+exp(-X1*t)/X1^2+(t-(t-tau).*heaviside(t-tau))/X1;
%D: convolution between exp(-X1*t) & t-(t-tau)H(t-tau)

A0_u=I1*(KuMr+4*MmMr/3)/((I0*MmMr-I1)^2)+KuMm/(I0*MmMr-I1);
B0_u=I1/(I0*MmMr-I1)-1/(MmMr-1);


u1=(-KlMr*delta*A0_u/tau)*D;
u0=(4/3)*KlMr*delta*(r_ratio^(2/3)-r_ratio^(-1/3))*(MmMr-1)*B0_u*D/(tau*(I0*MmMr-I1));


% 
% S2=-KuMr/(KuMr-4*MmMr/3);
% 
% 
% A0_u=KlMr*delta*KuMr/(4*(KuMr+KlMr*(r_ratio-1))+3*KlMr*KuMr*r_ratio); 
% A1_u=-KlMr*delta*((KuMr+4*MmMr/3)/(MmMr*I0-I1)+KuMr/I1);
% A0_Phi=1/I1;
% A1_Phi=((KuMr+4*MmMr/3)/(MmMr*I0-I1)+KuMr/I1)...
%     *(I1*MmMr*(MmMr-1)-MmMr*(I0*MmMr-I1))/(I1*MmMr*(KuMr-4*MmMr/3)+KuMr*MmMr*(I0*MmMr-I1));
% A2_Phi=(8/3)*MmMr*KuMm...
%     *(1/KuMr+(MmMr-1)/(KuMr-4*MmMr/3))/(I1*(KuMr-4*MmMr/3)+KuMr*(I0*MmMr-I1));
% 
% u1=A0_u+A1_u*exp(S1.*t);
% Phi=A0_Phi+A1_Phi*exp(S1.*t)+A2_Phi*exp(S2.*t);
% 
% u0=(4*KlMr*delta/3)*(r_ratio^(2/3)-r_ratio^(-1/3))*Phi+r_ratio^(2/3)*u1;
% Pl=-4*KlMr*KlMr*delta*(r_ratio-1)*Phi-3*KlMr*r_ratio*u1+KlMr*delta;
% 

