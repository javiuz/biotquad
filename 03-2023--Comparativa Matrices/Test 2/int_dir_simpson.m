function ints_g=int_dir_simpson(i,j,k,l,tt,comp)

global x y

ints_g=1/6*(sol_exactax(x(i,j),y(i,j),tt,comp)+...
    4*sol_exactax((x(i,j)+x(k,l))/2,(y(i,j)+y(k,l))/2,tt,comp)+...
    sol_exactax(x(k,l),y(k,l),tt,comp));

return
end 