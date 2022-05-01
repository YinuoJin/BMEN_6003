% Q3

%initialization
clear all;

%parameters
tspan=0:0.1:50; %time range

% oscillation setup
%{
a=30;
gx=200; 
gy=0.001;
gxy=0.001; 
m=0.01;
%}

% non-oscillation setup
a=50;
gx=0.5; 
gy=0.05;
gxy=0.1; 
m=1;

%% 3.3
y0 = [0,0];
[t,y] = ode23(@(t,y)posneg_fn(t,y,a,m,gx,gy,gxy), tspan, y0); %run ODE

figure(1);
plot(t, y(:,1), '-');
hold on;
xlabel('Time (t)');
ylabel('x');
%title('Oscillate');
title('Non-Oscillate');

%% 3.4
% phase plot
figure(2);
plot(y(:,1), y(:,2)); 
hold on;
xlabel('x');
ylabel('y');
title('Phase Portrait');


