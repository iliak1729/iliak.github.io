clear;
close all;
clc;
tic
%% Parameters
N = 301;
L = 1;
h = L/(N-1);
xset = 0:h:L;
yset = 0:h:L;

dt = .00002;
beta = 4.5;
Re = 5000;

%% Initialize
u_final = zeros(N);
v_final = zeros(N);
p_final = zeros(N);

% Two Sets of Staggered Variables
u = zeros(N+1,N);
v = zeros(N,N+1);
p = zeros(N+1)+1;
u(1,:) = 2;

u_new = zeros(N+1,N);
v_new = zeros(N,N+1);
p_new = zeros(N+1);
u_new(1,:) = 2;

%% Solve Equations Iteratively
error = 1;
iter = 0;
maxIter = Inf;
tol = 1e-8;
residuals = figure;

while (error > tol)

    % X-Momentum Interior
    for i = 2:N
        for j = 2:N-1
            uCx = ((.5*(u(i,j)+u(i,j+1)))^2-(.5*(u(i,j)+u(i,j-1)))^2)/h;
            uCy = (.25*(v(i-1,j)+v(i-1,j+1))*(u(i,j)+u(i-1,j))-.25*(v(i,j)+v(i,j+1))*(u(i,j)+u(i+1,j)))/h;
            uP = (p(i,j)-p(i,j+1))/h;
            uD = (u(i,j-1)-2*u(i,j)+u(i,j+1)+u(i-1,j)-2*u(i,j)+u(i+1,j))/(h*h*Re);

            dUdt = uP+uD-uCx-uCy;
            u_new(i,j) = u(i,j)+dt*dUdt;
        end
    end

    % X-Momentum Boundary Conditions
    u_new(1,:) = 2-u_new(2,:);
    u_new(N+1,:) = -u_new(N,:);
    u_new(2:N,1) = 0;
    u_new(2:N,N) = 0;

    % Y-Momentum Interior
    for i = 2:N-1
        for j = 2:N
            vCx = (.25*(u(i,j)+u(i+1,j))*(v(i,j)+v(i,j+1))-.25*(u(i,j-1)+u(i+1,j-1))*(v(i,j)+v(i,j-1)))/h;
            vCy = ((.5*(v(i,j)+v(i-1,j)))^2-(.5*(v(i,j)+v(i+1,j)))^2)/h;
            vP = (p(i+1,j)-p(i,j))/h;
            vD = (v(i,j-1)+v(i,j+1)+v(i+1,j)+v(i-1,j)-4*v(i,j))/(h*h*Re);

            dVdt = vP+vD-vCx-vCy;
            v_new(i,j) = v(i,j) + dt*dVdt;
        end
    end

    % Y-Momentum Boundary Conditions
    v_new(:,1) = - v_new(:,2);
    v_new(:,N+1) = - v_new(:,N);
    v_new(1,2:N) = 0;
    v_new(N,2:N) = 0;

    % Pressure Interior
    for i = 2:N
        for j = 2:N
            dUdx = (u(i,j)-u(i,j-1))/h;
            dVdy = (v(i-1,j)-v(i,j))/h;

            dPdt = -beta*(dUdx+dVdy);

            p_new(i,j) = p(i,j) + dt*dPdt;
        end
    end

    % Pressure Boundary Conditions
    p_new(1,:) = p_new(2,:);
    p_new(N+1,:) = p_new(N,:);
    p_new(:,N+1) = p_new(:,N);
    p_new(:,1) = p_new(:,2);

    %% Update Error and Show
    error = 0;
    for i = 2:N-1
        for j = 2:N-1
            dUdx = (u_new(i,j)-u_new(i,j-1))/h;
            dVdy = (v_new(i-1,j)-v_new(i,j))/h;

            error = error + abs(dUdx+dVdy);
        end
    end
    
    % Plot
    if(rem(iter,5000) == 0)
        figure(residuals);
        semilogy(iter,error,'-ko');
        hold on;
        xlabel("Iterations");
        ylabel("Error");
    end

    %% End Iteration
    u = u_new;
    v = v_new;
    p = p_new;
    iter = iter+1;
    % disp(iter);
end
%% Calc Final Vals
u_final = .5*(u(1:N,:)+u(2:N+1,:));
v_final = .5*(v(:,1:N)+v(:,2:N+1));
p_final = .25*(p(1:N,1:N)+p(2:N+1,1:N)+p(2:N+1,2:N+1)+p(1:N,2:N+1));
dVdx = (v(:,2:N+1)-v(:,1:N))/(h);
dUdy = (u(1:N,:)-u(2:N+1,:))/h;
vorticity = (dVdx-dUdy);
%% Visualize
x_dom = ((1:N)-1).*h;
y_dom = 1-((1:N)-1).*h;
[X,Y] = meshgrid(x_dom,y_dom);
figure;
hold on;
contourf(X,Y,u_final,21,'LineStyle','none');
xlabel("X");
ylabel("Y");
title("U");
colormap('jet');
colorbar;

figure;
hold on;
contourf(X,Y,v_final,21,'LineStyle','none');
xlabel("X");
ylabel("Y");
title("V");
colormap('jet');
colorbar;

figure;
hold on;
contourf(X,Y,p_final,21,'LineStyle','none');
xlabel("X");
ylabel("Y");
title("Pressure");
colormap('jet');
colorbar;

figure;
hold on;
xlabel("X");
ylabel("Y");
title("Vector Field");
quiver(X,Y,u_final,v_final,2.5);

figure;
hold on;
xlabel("X");
ylabel("Y");
title("Vorticity Field");
contourf(X,Y,vorticity,31,'LineStyle','none');
colorbar;
colormap('jet');

t = toc;

