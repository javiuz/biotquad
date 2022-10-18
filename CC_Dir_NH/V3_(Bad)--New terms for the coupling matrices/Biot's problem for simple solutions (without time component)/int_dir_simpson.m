function ints_g=int_dir_simpson(i,j,k,l,comp)

global x y

ints_g=1/6*(sol_exactax(x(i,j),y(i,j),comp)+...
    4*sol_exactax((x(i,j)+x(k,l))/2,(y(i,j)+y(k,l))/2,comp)+...
    sol_exactax(x(k,l),y(k,l),comp));

return
end 