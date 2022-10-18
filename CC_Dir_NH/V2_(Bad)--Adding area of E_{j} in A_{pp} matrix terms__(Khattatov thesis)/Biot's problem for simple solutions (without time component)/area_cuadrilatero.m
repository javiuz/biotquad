function area=area_cuadrilatero(x1,y1,x2,y2,x3,y3,x4,y4)

area=1/2*(abs(x2*y3+x3*y1+x1*y2-(x2*y1+x1*y3+x3*y2))+...
    abs(x4*y3+x3*y1+x1*y4-(x4*y1+x1*y3+x3*y4)));
return 
end