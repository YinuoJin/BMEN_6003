%positive/negative feedback function file
%source: class note

function dydt = pos_feedback_fn(t,y,a,g)

%y(1) is the first variable

dydt = zeros(1,1);
dydt(1) = (1+a*y(1)*y(1))/(1+y(1)*y(1))-g*y(1);



