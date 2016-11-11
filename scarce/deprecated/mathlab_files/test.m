clear;
dt = 0.01;
t=0:dt:10;
x0 = 0;
Q = 0;
v0 = 1;
x = x0;
x_vec = linspace(0, 0, size(t,2));
v_vec = linspace(0, 0, size(t,2));
dq_vec = linspace(0, 0, size(t,2));
q_vec = linspace(0, 0, size(t,2));
e_vec = linspace(0, 0, size(t,2));

index = 1;

for(i = t)
    v = v0+x*0.5*sin(0.1*x);
    dx = v * dt;
    x = x + dx;
    e = cos(0.2*x)+2;
    dq = e*dx;
    x_vec(index) = x;
    v_vec(index) = v;
    e_vec(index) = e;
    dq_vec(index) = dq;
    Q = Q + dq;
    q_vec(index) = Q;
    index = index + 1;
end;

plot(t,x_vec);
hold on;
plot(t,v_vec, 'COLOR', 'red');
plot(t,e_vec, 'COLOR', 'black');
plot(t,q_vec, 'COLOR', 'green');
xlabel('time [s]', 'FontWeight','bold');
ylabel('position/velocity', 'FontWeight','bold');
set(gcf, 'Color', [1 1 1]);
%ylim([1e-1 1e4]);
%hleg1 = legend('p-type FZ silicon','p-type oxigenated FZ silicon', 'n-type FZ silicon','n-type oxigenated FZ silicon', 'Location', 'northwest');
grid on;
set(gca, 'GridLineStyle', '-');
hold off;
    