clc;clear all;
%% Initialization
L = 1;
W = 1;
Nx = 10;
Ny = 10;
dx = L/Nx;
dy = W/Ny;
dt = 0.0001;
Re = 100;

u(1:Nx+2,1:Ny+2) = 0;
v(1:Nx+2,1:Ny+2) = 0;
psi(1:Nx+2,1:Ny+2) = 0;
omega(1:Nx+2,1:Ny+2) = 0;
omega_new = omega;
ap1 = -2*((1/(dx^2))+(1/(dy^2))); 

u_top = 0;
u_bottom = 1;
u_left = 0;
u_right = 0;

v_top = 0;
v_bottom = 0;
v_right = 0;
v_left = 0;

psi_top = 0;
psi_bottom = 0;
psi_left = 0;
psi_right = 0;
rms =0 ;

time_iteration_limit = 100;
loop_iteration_limit = 1;

%% Boundary Condition

%For u velocity
u(1,1:Ny+2) = 2*u_top - u(2,1:Ny+2);    
u(Nx+2,1:Ny+2) = 2*u_bottom - u(Nx+1,1:Ny+2);
u(1:Nx+2,1) = 2*u_left - u(1:Nx+2,2);
u(1:Nx+2,Ny+2) = 2*u_right - u(1:Nx+2,Ny+1);

%For v Velocity
v(1,1:Ny+2) = 2*v_top - v(2,1:Ny+2);    
v(Nx+2,1:Ny+2) = 2*v_bottom - v(Nx+1,1:Ny+2);
v(1:Nx+2,1) = 2*v_left - v(1:Nx+2,2);
v(1:Nx+2,Ny+2) = 2*v_right - v(1:Nx+2,Ny+1);

%For psi 
psi(1,1:Ny+2) = 2*psi_top - psi(2,1:Ny+2);  
psi(Nx+2,1:Ny+2) = 2*psi_bottom - psi(Nx+1,1:Ny+2);
psi(1:Nx+2,1) = 2*psi_left - psi(1:Nx+2,2);
psi(1:Nx+2,Ny+2) = 2*psi_right - psi(1:Nx+2,Ny+1);

%Updating the omega boundary Condition
omega_top = -8*((psi(2,1:Ny+2)-psi_top)/(dy^2)) - ((4*u_top)/(dy));
omega_bottom = -8*((psi(Nx+1,1:Ny+2)-psi_bottom)/(dy^2)) - ((4*u_bottom)/(dy));
omega_left =  -8*((psi(1:Nx+2,2)-psi_left)/(dx^2)) + ((4*u_left)/(dx));
omega_right =  -8*((psi(1:Nx+2,Ny+1)-psi_bottom)/(dx^2)) + ((4*u_right)/(dx));

%For omega
omega(1,1:Ny+2) = 2*omega_top - omega(2,1:Ny+2);    
omega(Nx+2,1:Ny+2) = 2*omega_bottom - omega(Nx+1,1:Ny+2);
omega(1:Nx+2,1) = 2*omega_left - omega(1:Nx+2,2);
omega(1:Nx+2,Ny+2) = 2*omega_right - omega(1:Nx+2,Ny+1);
omega
%% The Time stepping loop
for n = 1:time_iteration_limit
    
    for GSloop = 1:100000
        
        omega_new = omega;
        %For all the internal Cells
        for i = 2:Nx+1
            for j = 2:Ny+1

                %Flux Calculation
                Fe = 0.5*(u(i+1,j)+u(i,j))*dy;
                Fw = 0.5*(u(i,j)+u(i-1,j))*(-dy);
                Fn = 0.5*(v(i,j+1)+v(i,j))*dx;
                Fs = 0.5*(v(i,j)+v(i,j-1))*(-dx);

                %Simple Upwind 
                if Fe < 0
                    omega_e = omega(i-1,j);
                else
                    omega_e = omega(i,j);
                end

                if Fn < 0
                    omega_n = omega(i,j+1);
                else
                    omega_n = omega(i,j);
                end

                if Fw < 0
                    omega_w = omega(i+1,j);
                else
                    omega_w = omega(i,j);
                end

                if Fs < 0
                    omega_s = omega(i,j-1);
                else
                    omega_s = omega(i,j);
                end

                %Convective Term
                Convective_term = (Fe*omega_e) + (Fw*omega_w) + (Fn*omega_n) + (Fs*omega_s);

                %Diffusive Term
                Fde = -(1/Re)*((omega(i+1,j)-omega(i,j))/(dx))*(dy);
                Fdn = -(1/Re)*((omega(i,j+1)-omega(i,j))/(dy))*(dx);
                Fds = -(1/Re)*((omega(i,j-1)-omega(i,j))/(dy))*(dx);
                Fdw = -(1/Re)*((omega(i-1,j)-omega(i,j))/(dx))*(dy);
                Diffusive_term = Fde+Fdn+Fds+Fdw;

                %Updating Omega
                R = -(dx*dy)*(omega_new(i,j)-omega(i,j))/(dt)+(-Diffusive_term - Convective_term);
                ap = (dx*dy)*((1/dt)+(2/Re)*((1/(dx^2))+(1/(dy^2)))) + max([Fe,0])+max([Fw,0])+max([Fn,0])+max([Fs,0]);
                omega_new(i,j) = omega(i,j) + (R/ap);
                rms = rms + R^2;

                %% Updating BC
                omega_top = -8*((psi(2,1:Ny+2)-psi_top)/(dy^2)) - ((4*u_top)/(dy));
                omega_bottom = -8*((psi(Nx+1,1:Ny+2)-psi_bottom)/(dy^2)) - ((4*u_bottom)/(dy));
                omega_left =  -8*((psi(1:Nx+2,2)-psi_left)/(dx^2)) + ((4*v_left)/(dx));
                omega_right =  -8*((psi(1:Nx+2,Ny+1)-psi_right)/(dx^2)) + ((4*v_right)/(dx));

                omega_new(1,1:Ny+2) = omega_new(2,1:Ny+2)+ 2*(omega_top - omega_new(2,1:Ny+2));   
                omega_new(Nx+2,1:Ny+2) = omega_new(Nx+1,1:Ny+2) + 2*(omega_bottom - omega_new(Nx+1,1:Ny+2));
                omega_new(1:Nx+2,1) = omega_new(1:Nx+2,2) + 2*(omega_left - omega_new(1:Nx+2,2));
                omega_new(1:Nx+2,Ny+2) = omega_new(1:Nx+2,Ny+1) + 2*(omega_right - omega_new(1:Nx+2,Ny+1));


            end
        end
        
        rms = sqrt(rms/(Nx*Ny));
        RMS_For_omega = rms;
        RMS_For_omega
        if rms > 4e-2
            rms = 0;
            omega = omega_new;
        else
            fprintf("\n GS loop Convergence Reached!!!! \n")
            break;
        end 
    end
    
    Convective_term;
    Diffusive_term;
    %break
    
    %Solution for psi
    Rms = 1;
    while Rms > 1e-6
        
        Rms = 0;
        for i = 2:Nx+1
            for j = 2:Ny+1
                Residue = (-omega_new(i,j))-((psi(i+1,j)+psi(i-1,j)-2*psi(i,j))/(dx^2)) - ((psi(i,j-1)+psi(i,j+1)-2*psi(i,j))/(dy^2));
                psi(i,j) = psi(i,j) + (Residue./ap1);
                Rms = Rms + (Residue^2);
            end
        end
        
        psi(1,1:Ny+2) = 2*psi_top - psi(2,1:Ny+2);    
        psi(Nx+2,1:Ny+2) = 2*psi_bottom - psi(Nx+1,1:Ny+2);
        psi(1:Nx+2,1) = 2*psi_left - psi(1:Nx+2,2);
        psi(1:Nx+2,Ny+2) = 2*psi_right - psi(1:Nx+2,Ny+1);
        
        Rms = sqrt(Rms/(Nx*Ny));
        Rms_for_Psi = Rms;
        Rms_for_Psi
    end
    
    %u,v Calculation
    for i = 2:Nx+1
        for j = 2:Ny+1
            u(i,j) = (psi(i,j+1) - psi(i,j-1))/(2*dy);
            v(i,j) = -(psi(i+1,j) - psi(i-1,j))/(2*dx);   
        end
    end
    
    %Updating the BC    
    %For u velocity
    u(1,1:Ny+2) = 2*u_top - u(2,1:Ny+2);    
    u(Nx+2,1:Ny+2) = 2*u_bottom - u(Nx+1,1:Ny+2);
    u(1:Nx+2,1) = 2*u_left - u(1:Nx+2,2);
    u(1:Nx+2,Ny+2) = 2*u_right - u(1:Nx+2,Ny+1);

    %For v Velocity
    v(1,1:Ny+2) = 2*v_top - v(2,1:Ny+2);    
    v(Nx+2,1:Ny+2) = 2*v_bottom - v(Nx+1,1:Ny+2);
    v(1:Nx+2,1) = 2*v_left - v(1:Nx+2,2);
    v(1:Nx+2,Ny+2) = 2*v_right - v(1:Nx+2,Ny+1);
 
    %Steady state convergence Checking
    Rms = 0;
    for i = 2:Nx+1
        for j = 2:Ny+1
            Rms = Rms + ((omega_new(i,j) - omega(i,j))^2);
        end
    end
    Rms = sqrt(Rms/(Nx*Ny));
    Rms_k(n) = Rms/dt;

    
    omega_old = omega;
    omega = omega_new;
    
    Rms_time = Rms/dt
    
    if (Rms/dt) < 1e-3
        fprintf("Steady State Reached!!!!")
        break;
    end
    
end

%% Plotting
figure(1)
contourf(u,50)
colorbar
title("U Velocity profile")

figure(2)
contourf(v,50)
colorbar
title("V Velocity Profile")

figure(3)
contourf(psi,50)
colorbar
title("Psi Profile")

figure(4)
contourf(omega,70)
colorbar
title("Omega Profile")
%% CenterLine data
figure(5)
vertical_center_line_x = (Ny/2)+1;
u_velocity = u(2:Nx+1,vertical_center_line_x);
plot(u_velocity,(2:Ny+1)*dy)
xlabel("U Velocity")
ylabel("Y")
title("U Velocity vs y")

figure(6)
Horizontal_center_line_y = (Nx/2)+1;
v_velocity = v(Horizontal_center_line_y,2:Ny+1);
plot(2:Nx+1,v_velocity)
xlabel("X")
ylabel("V Velocity")
title("X vs V Velocity")