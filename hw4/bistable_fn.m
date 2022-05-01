function dydt = bistable_fn(t,y,a1,a2,b,g)
    dydt =[ 0; 0];
    x=y(1); y=y(2);
    
    dydt(1) = a1 / (1+y^b) - x;
    dydt(2) = a2 / (1+x^g) - y;

end