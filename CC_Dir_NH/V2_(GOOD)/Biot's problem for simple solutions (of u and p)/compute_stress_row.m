function [a1,b1,g1,a2,b2,g2,r,s]=compute_stress_row(sigma,row)

u11=sigma(2+row);
u12=sigma(row);
u21=sigma(4+row);
u22=sigma(6+row);
u31=sigma(8+row);
u32=sigma(10+row);
u41=sigma(14+row);
u42=sigma(12+row);

a1=-u11+u21+(-u22+u12+u32-u42)/2.;
b1=-u11+u41;
g1=u11;
a2=-u12+u22;
b2=-u12+u42+(-u21+u11-u41+u31)/2.;
g2=u12;
r=(u42-u32-u12+u22)/2.;
s=(u31-u21-u41+u11)/2.;

return
end