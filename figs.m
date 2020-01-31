t = linspace(-1,2.5,100);
[X,Y] = meshgrid(t,t);

Z = 100*(X.^2-Y).^2 + (X-1).^2;
%%
h = 0.005;
tspan = [0 20];
% y0 = [1.5;-1.2];%[.5;.5];%[3; 2];
% y0 =[-0.5436;0.6708];
%y0 = [0;0]
%% ODE45
f_R = @(x,y) 100*(x.^2-y).^2 + (x-1).^2;
gradf_R = @(t,y) -[400*y(1)*(y(1)^2-y(2))+2*(y(1)-1);-200*(y(1)^2-y(2))];
Mx = @(t,x) [400*x(1)^2+1, -200*x(1);-200*x(1), 100];
options = odeset('Mass',Mx);
figure;
hold on;
for i = 1:3
    y0 = 2.5*rand(1,2);
    plot3(y0(1),y0(2),f_R(y0(1),y0(2)),'k.','markersize',20)
    tspan = [0 50];
    
    [tout1,yout1] = ode45(gradf_R,tspan,y0);
    plot3(yout1(:,1),yout1(:,2),f_R(yout1(:,1),yout1(:,2)),'b','linewidth',2)
    plot3(yout1(end,1),yout1(end,2),f_R(yout1(end,1),yout1(end,2)),'b.','markersize',20)

    
    [tout2,yout2] = ode45(gradf_R,tspan,y0,options);
    plot3(yout2(:,1),yout2(:,2),f_R(yout2(:,1),yout2(:,2)),'r','linewidth',2);
    plot3(yout2(end,1),yout2(end,2),f_R(yout2(end,1),yout2(end,2)),'r.','markersize',20)
    
%     max(find(all(abs(yout1-[1;1])<1e-3),1),0)
    it_rg(i) = max(find(all(abs(yout1'-[1;1])<5e-3),1),0);
    it_ng(i) = max(find(all(abs(yout2'-[1;1])<5e-3),1),0);

end
contour3(X,Y,Z,500);
xlabel('x'); ylabel('y'); zlabel('z')
hold on; grid on
% axis([0 2.5 0 2.5 0 1000])
axis([-1 2.5 -1 2.5 0 1000])
plot3([1 1],[0 1],[0 0],'g--','linewidth',1.7)
plot3([0 1],[1 1],[0 0],'g--','linewidth',1.7)
set(gca,'fontname','times','fontsize',12,'fontweight','bold');
%% Forward EULER
h = [0.001,0.002,0.005,0.01,0.05];

tspan = [0 20];
y0 = [1.5;.5];%[3; 2];
gradf_R = @(t,y) -[400*y(1)*(y(1)^2-y(2))+2*(y(1)-1);-200*(y(1)^2-y(2))];
natGrad = @(t,x) Mx(t,x)\gradf_R(t,x);

figure;
plot3(y0(1),y0(2),f_R(y0(1),y0(2)),'k.','markersize',20)
hold on;
for i = 1:2%length(h)
    
    [tout1,yout1] = forwardEulerSys(gradf_R,tspan,y0,h(i)); yout1 = yout1';
    
    plot3(yout1(:,1),yout1(:,2),f_R(yout1(:,1),yout1(:,2)),'k','linewidth',2.5)
    plot3(yout1(end,1),yout1(end,2),f_R(yout1(end,1),yout1(end,2)),'b.','markersize',20)
    
    [tout2,yout2] = forwardEulerSys(natGrad,tspan,y0,h(i)); yout2 = yout2';
    plot3(yout2(:,1),yout2(:,2),f_R(yout2(:,1),yout2(:,2)),'r','linewidth',2.5);
    plot3(yout2(end,1),yout2(end,2),f_R(yout2(end,1),yout2(end,2)),'r.','markersize',20)
end
contour3(X,Y,Z,500);
xlabel('x'); ylabel('y'); zlabel('z')
hold on; grid on
axis([-1 2.5 -1 2.5 0 1000])
plot3([1 1],[-1 1],[0 0],'g--','linewidth',1.7)
plot3([-1 1],[1 1],[0 0],'g--','linewidth',1.7)
set(gca,'fontname','times','fontsize',13,'fontweight','bold');
%%
axis([-1 2.5 -1 2.5 0 1000])