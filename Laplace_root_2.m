function y=Laplace_root_2(x,r_o,h0,f1)
%This function plots equation (A.57), which will be used to identify the
%postinjection time scale
% the input: r_o=r_o/R_o; h0 and f1 are defined in Appendix
%This code was used to compute results shown in Figure A.2 (b)
y=sin(sqrt(x)*(1-r_o)).*((f1/r_o-h0*r_o)*x-h0)-sqrt(x).*(f1*(x)/r_o+h0*r_o-h0).*cos(sqrt(x)*(1-r_o));
%f=sin(sqrt(x)*(1-ro)).*((ho*ro-f1/ro).*x+ho)-cos(sqrt(x)*(1-ro)).*sqrt(x).*(ho*(1-ro)-f1*x/ro);
end
%[300,80,50,30,30,30,20,20,20]