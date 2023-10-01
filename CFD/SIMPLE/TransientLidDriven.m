clear;
close all;
clc;
%% Input Parameters
Nx = 51;
Ny = 51;
Re = 100;
cd 'C:\Users\ik465\Desktop\School\Personal\SIMPLE\Videos'
dt = 1000000;
N = 1; % Must be divisible evenly by frames
frames = 1;
tol = 1e-6;
fps = 30;
name = "Transient Lid Driven Flow - Re "+num2str(Re)+' - '+num2str(Nx) + 'x'+num2str(Ny) +' - ' + num2str(N)+"x"+num2str(dt)+ 's Sim' + num2str(frames) + ' Frames '+ num2str(fps) + 'fps'
window = 10000; % For Residual Viewing
% Under-relaxation
alpha = .8;
timeLength = 0;
%% Derived Parameters
xset = linspace(0,1,Nx);
dx = xset(2) - xset(1);
yset = linspace(0,1,Ny);
dy = yset(2) - yset(1);
[X,Y] = meshgrid(xset,1-yset);


%% Allocate Space for Variables
% Saved Colocated Values
u_saved = zeros(Ny,Nx,frames+1);
v_saved = zeros(Ny,Nx,frames+1);
p_saved = zeros(Ny,Nx,frames+1);

% Staggered Variables
u = zeros(Ny+1,Nx);
u(1,:) = 2;
u_star = u;
u_m = u;
d_e(Ny+1,Nx)=0;

v = zeros(Ny,Nx+1);
v_star = zeros(Ny,Nx+1);
v_m = zeros(Ny,Nx+1);
d_n=zeros(Ny,Nx+1);

p = zeros(Ny+1,Nx+1);
p_m = zeros(Ny+1,Nx+1);
pc = zeros(Ny+1,Nx+1);
b = zeros(Ny+1,Nx+1);

% Initial Conditions
u0 = u;
v0 = v;
p0 = p;

u_prev = u0;
v_prev = v0;
p_prev = p0;
% Save Initial Condition Frame
for i = 1:Ny
    for j = 1:Nx
        u_saved(i,j,1) = 0.5*(u0(i,j) + u0(i+1,j));
        v_saved(i,j,1) = 0.5*(v0(i,j) + v0(i,j+1));
        p_saved(i,j,1) = 0.25*(p0(i,j) + p0(i,j+1) + p0(i+1,j) + p0(i+1,j+1));
    end
end
%% Solve Governing Equations
iter = 0;
count = 1;
for t = 1:N
    graphError = 1;
    error = 1;
    in = 1;
    tic;
    while error > tol || in < 3

        % X Momentum
        for i = 2:Ny
            for j = 2:Nx-1
                u_e = (u(i,j) + u(i,j+1))/2;
                u_w = (u(i,j) + u(i,j-1))/2;
                v_n = (v(i-1,j) + v(i-1,j+1))/2;
                v_s = (v(i,j) + v(i,j+1))/2;

                A_E = u_e*dy/2-dy/(Re*dx);
                A_W = -u_w*dy/2-dy/(Re*dx);
                A_N = v_n*dx/2-dx/(Re*dy);
                A_S = -v_s*dx/2-dx/(Re*dy);

                A_P = (u_e*dy-u_w*dy+v_n*dx-v_s*dx)/2+2*dy/(Re*dx)+2*dx/(Re*dy)+dx*dy/dt;
                d_e(i,j) = -dy/A_P;
                
                Q_unsteady = u_prev(i,j)*dx*dy/dt;
                u_star(i,j) = (Q_unsteady-A_E*u(i,j+1)-A_W*u(i,j-1)-A_N*u(i-1,j)-A_S*u(i+1,j))/A_P - d_e(i,j) * (p(i,j) - p(i,j+1));
            end
        end
        % X Momentum Boundary
        u_star(1,:) = 2 - u_star(2,:); % Top BC
        u_star(Ny + 1,:) = -u_star(Ny,:); % Bottom BC
        u_star(2:Ny,1) = 0; % Left BC
        u_star(2:Ny,Nx) = 0; % Right BC

        % Y Momentum
        for i = 2:Ny-1
            for j = 2:Nx
                u_w = (u(i,j-1) + u(i+1,j-1))/2;
                u_e = (u(i,j) + u(i+1,j))/2;
                v_n = (v(i,j) + v(i-1,j))/2;
                v_s = (v(i,j) + v(i+1,j))/2;

                A_E = u_e*dy/2-dy/(Re*dx);
                A_W = -u_w*dy/2-dy/(Re*dx);
                A_N = v_n*dx/2-dx/(Re*dy);
                A_S = -v_s*dx/2-dx/(Re*dy);

                A_P = (u_e*dy-u_w*dy+v_n*dx-v_s*dx)/2+2*dy/(Re*dx)+2*dx/(Re*dy)+dx*dy/dt;
                d_n(i,j) = -dx/A_P;

                Q_unsteady = v_prev(i,j)*dx*dy/dt;
                v_star(i,j) = (Q_unsteady-A_E*v(i,j+1)-A_W*v(i,j-1)-A_N*v(i-1,j)-A_S*v(i+1,j))/A_P - d_n(i,j) * (p(i+1,j) - p(i,j));
            end
        end

        % Y Momentum Boundary
        v_star(:,1) = -v_star(:,2);
        v_star(:,Nx + 1) = -v_star(:,Nx);
        v_star(1,2:Ny) = 0;
        v_star(Ny,2:Nx) = 0;
        
        % Pressure Correction
        pc(1:Ny+1,1:Nx+1)=0;
        for i = 2:Ny
            for j = 2:Nx
                A_E = -d_e(i,j)*dy;
                A_W = -d_e(i,j-1)*dy;
                A_S = -d_n(i,j)*dx;
                A_N = -d_n(i-1,j)*dx;
                A_P = A_E+A_W+A_S+A_N;

                b(i,j) = -(u_star(i,j) - u_star(i,j-1))*dy + (v_star(i,j) - v_star(i-1,j))*dx;
            
                pc(i,j) = (A_E*pc(i,j+1) + A_W*pc(i,j-1) + A_N*pc(i-1,j) + A_S*pc(i+1,j) + b(i,j))/A_P;
                pc(i,j) = b(i,j)/A_P;
            end
        end

        % Apply Pressure Correction
        for i = 2:Ny
            for j = 2:Nx
                p_m(i,j) = p(i,j) + alpha*pc(i,j);
            end
        end
        
        % Pressure Boundary
        p_m(1,:) = p_m(2,:);
        p_m(Ny + 1,:) = p_m(Ny,:);
        p_m(:,1) = p_m(:,2);
        p_m(:,Nx + 1) = p_m(:,Nx);

        % Correcting the velocities

        % X Correction
        for i = 2:Ny
            for j = 2:Nx - 1
                u_m(i,j) = u_star(i,j) + alpha*d_e(i,j)*(pc(i,j+1) - pc(i,j));
            end
        end
        
        % X Boundary
        u_m(1,:) = 2 - u_m(2,:);
        u_m(Ny + 1,:) = -u_m(Ny,:);
        u_m(2:Ny,1) = 0;
        u_m(2:Ny,Nx) = 0;
        
        % Y Correction
        for i = 2:Ny - 1
            for j = 2:Nx
                v_m(i,j) = v_star(i,j) + alpha*d_n(i,j)*(pc(i,j) - pc(i+1,j));
            end
        end

        % Y Boundary
        v_m(:,1) = -v_m(:,2);
        v_m(:,Nx + 1) = -v_m(:,Nx);
        v_m(1,2:Nx) = 0;
        v_m(Ny,2:Nx) = 0;

        % Calculate Error
        error = 0;
        for i = 2:Ny
            for j = 2:Nx
                error = error + abs(b(i,j));
            end
        end

        if(in == 1)
            graphError = error;

            % disp(error);
        end
        in = in +1;
        % Update Residual Graph
%         if(rem(iter, 100) == 0 && iter > 0)
%            prevError = graphError;
%            prevIter = iter - 100;
%            graphError = error;
%            figure(1);
%            %semilogy([prevIter iter], [prevError error], '-ko') % 304.6752
%            semilogy(iter,error,'ko'); % 286.467
%            hold on
%            xlabel('Iterations')
%            ylabel('Residual Error')
%            if(iter > window)
%                 xlim([iter-window, iter])
%            end
%         end
        u = u_m;
        v = v_m;
        p = p_m;
        iter = iter + 1;
    end
    time = toc;
    % Once Converged, save variables for frame if necessary
    if(mod(t,round(N/frames)) == 0)
        count = count + 1;
        for i = 1:Ny
            for j = 1:Nx
                u_saved(i,j,count) = 0.5*(u(i,j) + u(i+1,j));
                v_saved(i,j,count) = 0.5*(v(i,j) + v(i,j+1));
                p_saved(i,j,count) = 0.25*(p(i,j) + p(i,j+1) + p(i+1,j) + p(i+1,j+1));
            end
        end
    end
    
    % Update Previous Timestep information
    u_prev = u_m;
    v_prev = v_m;
    p_prev = p_m;
    u = u_m;
    v = v_m;
    p = p_m;

    disp("=============================================");
    disp(name);
    disp(num2str(t)+"/" + num2str(N));
    disp("Timestep computaiton time: " + num2str(time) +" Seconds");
    remIters = N-t;
    timeLeft = remIters*time;
    secLeft = mod(timeLeft,60);
    timeLeft = timeLeft-secLeft;
    timeLeft = timeLeft/60;
    minLeft = mod(timeLeft,60);
    timeLeft = timeLeft - minLeft;
    timeLeft = timeLeft/60;
    disp("Approx. Time Left: "+num2str(timeLeft)+"h"+num2str(minLeft)+"min"+num2str(secLeft)+"sec");
    disp("Converged Error: "+num2str(error));
    disp("=============================================");
    timeLength = timeLength + time;
    error = 1;
end

%% Visualize Data
 % create the video writer with 1 fps
 % set the seconds per image
 secsPerImage = [5 10 15];
 % open the video writer

 writerObj = VideoWriter(name,"MPEG-4");
 writerObj.FrameRate = 60;
 open(writerObj);

figure(2);
lengthPrev = -1;
for i = 1:frames+1
    t = dt*(i-1)*N/frames;
    contourf(X,Y,u_saved(:,:,i), 21, 'LineStyle', 'none')
    tStr = num2str(t,'%05.2f');
    title("U Velocity Contour - t = "+tStr+"s");
    colorbar
    colormap('jet')
    xlabel('x')
    ylabel('y')
    % Write Video
    F = getframe(gcf);
    writeVideo(writerObj, F);
    lengthPrev = length(tStr);
end

close(writerObj);
disp("EXPORT COMPLETE - NAME:");
disp(name);

