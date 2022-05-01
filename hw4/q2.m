% Q2
% Code modified from class note

%initialization
clear all;

%parameters

a=50; %strength of pos feedback
x = 0:0.01:10;
g_vals=[10, 13, 16, 19, 22, 25]; %decay rate
dx = zeros(length(g_vals), length(x));
r = zeros(6, 3);

%% 2.2 
%figure(1);
for i=1:length(g_vals)
    g = g_vals(i);
    dx(i,:) = (1+a.*x.^2) ./ (1+x.^2) - g*x;
    plot(x, dx(i,:), '.-'); hold on; %plot
    
    % get real parts of fixed points
    r(i,:) = sort(real(roots([g, -a, g, -1])));

end
plot(x, 0, '.-');
legend('g=10', 'g=13', 'g=16', 'g=19', 'g=22', 'g=25');
xlabel('x');
ylabel('dx');


%% 2.3
% print out roots

%for i=1:size(r, 1)
%    disp(r(i,:);
%end

%% 2.4
mat_is_stable = zeros(6, 3);
for i=1:size(r, 1)
    for j=1:size(r, 2)
        x_fix = r(i,j);
        mat_is_stable(i,j) = 2*x_fix*(a-1) / (1+x_fix^2)^2 - g;
    end
end


