function dydt = posneg_fn(t,y,a,m,gx,gy,gxy)
    dydt =[ 0; 0];
    x=y(1); y=y(2);
    
    dydt(1) = (1+a*x^2) / (1+x^2) - gx*x - gxy*xy;
    dydt(2) = m*(1+a*x^2) / (1+x^2) - gy*y;

end


