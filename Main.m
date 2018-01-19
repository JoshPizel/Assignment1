%clear all

clearvars
clearvars -GLOBAL
close all

global C
global X Y

    C.q_0 = 1.60217653e-19;             % electron charge
    C.hb = 1.054571596e-34;             % Dirac constant
    C.h = C.hb * 2 * pi;                % Planck constant
    C.m_0 = 9.10938215e-31;             % electron mass
    C.kb = 1.3806504e-23;               % Boltzmann constant
    C.eps_0 = 8.854187817e-12;          % vacuum permittivity
    C.mu_0 = 1.2566370614e-6;           % vacuum permeability
    C.c = 299792458;                    % speed of light
    C.g = 9.80665;                      % metres (32.1740 ft) per s²
    
mn=0.26*C.m_0; %electron mass
Temp = 300; %Given in kelvin
rTime=1000; %run time in timesteps

%thermal velocity
Vth = sqrt(C.kb*Temp/mn);

%establish inital electron positions
%working area 200nm x 100nm
workX=200*10^-9;
workY=100*10^-9;

size=10;

X= rand(3,size);
Y= rand(3,size);

%positions
X(1,:)= X(1,:)*workX;
Y(1,:)= Y(1,:)*workY;

%initial direction of each particle
X(3,:) = X(3,:)*pi;
Y(3,:) = Y(3,:)*pi;

%velocity of each particle
for i =1:length(X)
    
    if(X(3,i)<pi/2)
        X(2,i) = Vth*X(2,i)*cos(X(3,i));
    else
        X(2,i) = Vth*X(2,i)*cos(X(3,i)-pi/2)*(-1);
    end
    
    if(Y(3,i)<pi/2)
        Y(2,i) = Vth*Y(2,i)*cos(Y(3,i));
    else
        Y(2,i) = Vth*Y(2,i)*cos(Y(3,i)-pi/2)*(-1);
    end
    
end

%set timestep of function
spacStep = 0.01*workY;
dt = 1/Vth*spacStep;
steps = rTime/dt;

%variable change
%setup mapping


%main function
for i = 0:dt:steps
    
    
end


figure(1)
scatter(X(1,:),Y(1,:))
xlim([0 workX])
ylim([0 workY])