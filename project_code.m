f_R = @(x) 100*(x(1)^2-x(2))^2 + (x(1)-1)^2;
%%
t = linspace(-2.5,2.5,100);
[X,Y] = meshgrid(t,t);

Z = 100*(X.^2-Y).^2 + (X-1).^2;

% figure;
% surf(X,Y,Z);

% 
figure;
contour3(X,Y,Z,500);

axis([-2.5 2.5 -2.5 2.5 0 1000])
%%
h = 0.005;
tspan = [0 20];
% y0 = [1.5;-1.2];%[.5;.5];%[3; 2];
% y0 =[-0.5436;0.6708];
%y0 = [0;0]
%% ODE45
y0 = [0 2];
tspan = [0 30];
%y0 = [2; 2];
gradf_R = @(t,y) -[400*y(1)*(y(1)^2-y(2))+2*(y(1)-1);-200*(y(1)^2-y(2))];
[tout1,yout1] = ode45(gradf_R,tspan,y0);

f_R = @(x,y) 100*(x.^2-y).^2 + (x-1).^2;
figure;
% surf(X,Y,Z);
hold on;
plot3(yout1(:,1),yout1(:,2),f_R(yout1(:,1),yout1(:,2)),'linewidth',2)
% axis([-2 2 -2 2 0 1000])
plot3(y0(1),y0(2),f_R(y0(1),y0(2)),'k.','markersize',20)
plot3(yout1(end,1),yout1(end,2),f_R(yout1(end,1),yout1(end,2)),'b.','markersize',20)

Mx = @(t,x) [400*x(1)^2+1, -200*x(1);-200*x(1), 100];
options = odeset('Mass',Mx);
[tout2,yout2] = ode45(gradf_R,tspan,y0,options);

plot3(yout2(:,1),yout2(:,2),f_R(yout2(:,1),yout2(:,2)),'r','linewidth',2);
plot3(yout2(end,1),yout2(end,2),f_R(yout2(end,1),yout2(end,2)),'r.','markersize',20)

contour3(X,Y,Z,500);
xlabel('x'); ylabel('y')
%% Forward EULER
h = 0.002;
tspan = [0 20];
y0 = [.5;.5];%[3; 2];
gradf_R = @(t,y) -[400*y(1)*(y(1)^2-y(2))+2*(y(1)-1);-200*(y(1)^2-y(2))];
[tout1,yout1] = forwardEulerSys(gradf_R,tspan,y0,h); yout1 = yout1';
figure;
%surf(X,Y,Z);
hold on;
plot3(yout1(:,1),yout1(:,2),f_R(yout1(:,1),yout1(:,2)),'linewidth',2)
% axis([-2 2 -2 2 0 1000])
plot3(y0(1),y0(2),f_R(y0(1),y0(2)),'k.','markersize',20)
plot3(yout1(end,1),yout1(end,2),f_R(yout1(end,1),yout1(end,2)),'b.','markersize',20)


natGrad = @(t,x) Mx(t,x)\gradf_R(t,x);
[tout2,yout2] = forwardEulerSys(natGrad,tspan,y0,h); yout2 = yout2';
plot3(yout2(:,1),yout2(:,2),f_R(yout2(:,1),yout2(:,2)),'r','linewidth',2);
plot3(yout2(end,1),yout2(end,2),f_R(yout2(end,1),yout2(end,2)),'r.','markersize',20)

xlabel('x'); ylabel('y')


contour3(X,Y,Z,500);
%% Extended Rosenbrock
fxN = @(x) 100*(x(1).^2-x(2)).^2 + 2*(x(1)-1).^2 + ...
            100*(x(2).^2-x(3)).^2 + 2*(x(2)-1).^2 + ...
            100*(x(3).^2-x(4)).^2 + 2*(x(3)-1).^2 ;
        
        
grad_fN = @(t,x) - [400*x(1)*(x(1)^2-x(2))+2*(x(1)-1);...
                    -200*(x(1)^2-x(2))+400*x(2)*(x(2)^2-x(3))+2*(x(2)-1);...
                    -200*(x(2)^2-x(3))+400*x(3)*(x(3)^2-x(4))+2*(x(3)-1);...
                    -200*(x(3)^2-x(4))];
                
h = 0.001;
tspan = [0 10];
y0 = [-1.5;1.2;-.5;-.1];%[.5;.5];%[3; 2];
y0 =[-0.5436;0.6708;0.0139;0.9245];
[tout1,yout1] = ode45(grad_fN,tspan,y0);

figure;
subplot(1,2,1)
plot(tout1,yout1(:,1));
hold on;
plot(tout1,yout1(:,2));
plot(tout1,yout1(:,3));
plot(tout1,yout1(:,4));

[tout1,yout1] = forwardEulerSys(grad_fN,tspan,y0,h);yout1 = yout1';

subplot(1,2,2)
plot(tout1,yout1(:,1));
hold on;
plot(tout1,yout1(:,2));
plot(tout1,yout1(:,3));
plot(tout1,yout1(:,4));
        
%% 
% y0 = 2*rand(4,1)-1;
y0 =[-0.5436;0.6708;0.0139;0.9245];
tspan = [0 20];

Theta = @(x) [20*x(1) -10 0 0;...
              1 0 0 0;...
              0 20*x(2) -10 0;...
              0 1 0 0;...
              0 0 20*x(3) -10;...
              0 0 1 0];
        
Mx = @(t,x) Theta(x)'*Theta(x);
options = odeset('Mass',Mx);
[tout2,yout2] = ode45(grad_fN,tspan,y0,options);
figure;
subplot(1,2,1)
plot(tout2,yout2(:,1));
hold on;
plot(tout2,yout2(:,2));
plot(tout2,yout2(:,3));
plot(tout2,yout2(:,4));     
        
h = 0.001;

natGrad = @(t,x) Mx(t,x)\grad_fN(t,x);
[tout2,yout2] = forwardEulerSys(natGrad,tspan,y0,h); yout2 = yout2';
subplot(1,2,2)
plot(tout2,yout2(:,1));
hold on;
plot(tout2,yout2(:,2));
plot(tout2,yout2(:,3));
plot(tout2,yout2(:,4));     
