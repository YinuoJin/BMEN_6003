% Q4

%initialization
clear all;

%parameters
tspan=0:0.01:50; %time range
a1 = 10;
a2 = 10;
b = 0.6;
g = 0.6;
y0 = [1,0];

[t,y] = ode23(@(t,y)bistable_fn(t,y,a1,a2,b,g), tspan, y0); %run ODE

figure(1);
plot(t, y(:,1), '-');
xlabel('Time (t)');
ylabel('x');
title('X - Bistable toggle switch');

figure(2);
plot(t, y(:,2), '-');
xlabel('Time (t)');
ylabel('y');
title('Y - Bistable toggle switch');

