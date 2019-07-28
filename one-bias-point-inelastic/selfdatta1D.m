function EEE = selfdatta1D(E,U,t0)

x = (E-U)/(2*t0);

if (0<=x)&&(x<=2)
    EEE = t0*(x-1)-1i*t0*sqrt(2*x-x^2);
elseif x<= 0
    EEE = t0*(x-1) + t0*sqrt(x^2 -2*x);
else x>=2
     EEE = t0*(x-1)-t0*sqrt(x^2 -2*x);
    
end  
    
