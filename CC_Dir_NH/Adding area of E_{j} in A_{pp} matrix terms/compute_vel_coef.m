function [a1,b1,g1,a2,b2,g2,r,s]=compute_vel_coef(vel)

u11=vel(2);
u12=vel(1);
u21=vel(3);
u22=vel(4);
u31=vel(5);
u32=vel(6);
u41=vel(8);
u42=vel(7);

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