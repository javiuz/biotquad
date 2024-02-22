function sig=exact_sol_sigma(xx,yy,tt)

global alpha mu lambda

s11=-alpha*tt*xx + 2*(lambda + mu); 
s12=0;
s22=-alpha*tt*xx + 2*(lambda + mu);   

sig=[s11 s12;s12 s22];

return
end
