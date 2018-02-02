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
rTime=10000; %run time in timesteps
MTBC = 0.2e-12;

%thermal velocity
Vth = sqrt(2*C.kb*Temp/mn);

%establish inital electron positions
%working area 200nm x 100nm
workX=200*10^-9;
workY=100*10^-9;

size=1000;
displaySize=10;

X= rand(2,size);
Y= rand(2,size);

%positions
X(1,:)= X(1,:)*workX;
Y(1,:)= Y(1,:)*workY;

colour = rand(1,displaySize);
%initial direction of each particle
angle(1,:) = X(2,:)*2*pi;

% for normal distribution of velocity
sigma =sqrt(C.kb*Temp/mn)/4;
mu = Vth;
MBdist = makedist('Normal',mu,sigma);
Xvel = zeros(1,size);
Yvel = zeros(1,size);
for i=1:1:size
    vel = random(MBdist);
    Xvel(1,i) = vel*cos(angle(1,i));
    Yvel(1,i) = vel*sin(angle(1,i));
    
    %check = velocity<0;
    %angle(check) = angle(check)+pi;
end

%Xvel(1,:) = Vth*cos(angle(1,:));
%Yvel(1,:) = Vth*sin(angle(1,:));

%hist(velocity)
%set timestep of function
spacStep = 0.01*workY;
dt = spacStep/Vth;
steps = 1000;

%variable change
%setup mapping
Xpos = X(1,:);
Ypos = Y(1,:);

Xvel(1,:) = Xvel(1,:)*dt;
Yvel(1,:) = Yvel(1,:)*dt;

%percent scatter
Pscat=1-exp(-(dt/MTBC));

figure(1)
%main function
for i = 1:1:steps
    
    %scattering  
    scattered=rand(1,size);
    scatterCheck = scattered<=Pscat;
    angle(scatterCheck) = rand*2*pi;
    velocity = random(MBdist,1,size);
    Xvel(scatterCheck) = velocity(scatterCheck).*cos(angle(scatterCheck))*dt;
    Yvel(scatterCheck) = velocity(scatterCheck).*sin(angle(scatterCheck))*dt;
    
    %position advance
    %logical indexing
    checkXright = Xpos +Xvel>2e-7;%right side period
    Xpos(checkXright) = Xpos(checkXright)+Xvel(checkXright)-workX;
    checkXleft = Xpos +Xvel<0;%left side period
    Xpos(checkXleft) = Xpos(checkXleft) +Xvel(checkXleft)+workX;
    
    %bottle neck
    checkXbox = (Xpos + Xvel)>0.8e-7 && (Xpos + Xvel)<1.2e-7;
    checkYBoxbottom =(Ypos + Yvel)<0.4e-7;
    checkBoxBottom = checkYBoxbottom && checkXbox;
    
    Xvel(checkBoxBottom) = Xvel(checkBoxBottom).*(-1);
    
    leftover = ~(checkXright | checkXleft);
    
    Xpos(leftover) = Xpos(leftover) +Xvel(leftover);
    
    checkY = (Ypos+Yvel>1e-7 | Ypos+Yvel<0);%reflection y's
    Yvel(checkY) = Yvel(checkY).*(-1);
    Ypos(1,:) = Ypos(1,:)+Yvel(1,:);
    
    %temperature calculations
    Ysum = sum((Yvel/dt).^2);
    Xsum = sum((Xvel/dt).^2);
    calcTemp = mn*((Ysum)+(Xsum))/(2*C.kb);
    averageTemp = calcTemp/size;
    
    %plotting here
    prevX(i,:) =Xpos(1,:);
    prevY(i,:) =Ypos(1,:);
    for j = 1:1:displaySize
        plot(prevX(:,j),prevY(:,j),'color',[colour(1,j) 0 j/displaySize])
        xlim([0 workX])
        ylim([0 workY])
        legend(['Temperature:' num2str(averageTemp)])
        drawnow
        hold on
    end
    
    
end