function [r_o,h0,f1,h1,M,tensile,t,xx,E1]=Transient_analytical_solver(t_range,dt,alpha,phi_o,r_ratio,KsMr,KfMr,KlMr,MmMr,delta)
% This is the code that solves for the evolution of total transported mass
% following a sudden injection, when only the slowest decay term in the
% analytical solution is included. The result corresponds to the blue curve
% in Figure A2 (a). 
KsMm=KsMr/MmMr;
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

r_o=r_ratio^(-1/3);

fun=@(x) Laplace_root_2(1/x,r_o,h0,f1);
%!! Only if time scale is smaller than 0.5\tau_diff
ts=fzero(fun,0.5); %this will give the largest value of 1/x, for root x
xx=1/ts;

E1=f1^2*(1-r_o)*xx^2/(2*r_o^2)...
    -(2*h0*f1-3*h0*f1*r_o+f1^2+h0*r_o^3*(r_o-1))*xx/(2*r_o)...
    +h0^2*(1-r_o^3)/2+1.5*h0*f1;
t=t_range(1):dt:t_range(2);
M=3*r_ratio*h1*((1-1/r_ratio)/(3*f1+h0*(1-1/r_ratio))-f1*(xx+1)*exp(-xx*t)/E1);
tensile=2*(-g1*M/(g2*3*r_ratio)-f4/g2);
