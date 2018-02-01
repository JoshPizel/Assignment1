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
<<<<<<< HEAD:Main.m
rTime=1000; %run time in timesteps
betCollisions =0.2e-9; %mean time between collisions
=======
rTime=10000; %run time in timesteps
>>>>>>> 4c268e812079f8c8927e653dda8e3e1d050bd937:assignment1/Main.m

%thermal velocity
Vth = sqrt(2*C.kb*Temp/mn);

%establish inital electron positions
%working area 200nm x 100nm
workX=200*10^-9;
workY=100*10^-9;

size=1000;
displaySize=10;

%set timestep of function
spacStep = 0.01*workY;
dt = spacStep/Vth;
steps = 1000;

%start probability density function for one dimension
sigma = sqrt(mn/(2*pi*C.kb*Temp));
MBdist = makedist('Normal','mu', Vth, 'sigma', sigma);
valueBin = linspace(Vth-3*sigma, Vth+3*sigma,100);

valueDistro = pdf(MBdist,valueBin);

histogram(valueDistro,100)

X= rand(2,size);
Y= rand(2,size);

%positions
X(1,:)= X(1,:)*workX;
Y(1,:)= Y(1,:)*workY;

colour = rand(1,displaySize);
%initial direction of each particle
angle(1,:) = X(2,:)*2*pi;

<<<<<<< HEAD:Main.m
%velocity of each particle
X(2,:) = Vth*cos(angle(1,:));
Y(2,:) = Vth*sin(angle(1,:));
=======
%for normal distribution of velocity
sigma =sqrt(C.kb*Temp/mn);
mu = Vth;
MBdist = makedist('Normal',mu,sigma);
velocity = zeros(1,size);
Xvel = zeros(1,size);
Yvel = zeros(1,size);
for i=1:1:size
    velocity(1,i) = random(MBdist);
    Xvel(1,i) = velocity(1,i)*cos(angle(1,i));
    Yvel(1,i) = velocity(1,i)*sin(angle(1,i));
end
%set timestep of function
hist(velocity)
spacStep = 0.01*workY;
dt = spacStep/Vth;
steps = 1000;
>>>>>>> 4c268e812079f8c8927e653dda8e3e1d050bd937:assignment1/Main.m

%variable change
%setup mapping
Xpos = X(1,:);
Ypos = Y(1,:);

Xvel(1,:) = Xvel(1,:)*dt;
Yvel(1,:) = Yvel(1,:)*dt;

MFP = betCollisions*Vth;

figure(1)
%main function
for i = 1:1:steps
    %position advance
    %logical indexing
    checkXright = Xpos +Xvel>2e-7;
    Xpos(checkXright) = Xpos(checkXright)+Xvel(checkXright)-workX;
    checkXleft = Xpos +Xvel<0;
    Xpos(checkXleft) = Xpos(checkXleft) +Xvel(checkXleft)+workX;
    
    leftover = ~(Xpos +Xvel>2e-7 | Xpos +Xvel<0);
    
    Xpos(leftover) = Xpos(leftover) +Xvel(leftover);
    
    checkY = (Ypos+Yvel>1e-7 | Ypos+Yvel<0);
    Yvel(checkY) = Yvel(checkY).*(-1);
    Ypos(1,:) = Ypos(1,:)+Yvel(1,:);
    
    Ysum = sum(Yvel);
    Xsum = sum(Xvel);
    calcTemp = mn*((Ysum/dt)^2+(Xsum/dt)^2)/(2*C.kb);
    averageTemp = calcTemp/size;
    
    prevX(i,:) =Xpos(1,:);
    prevY(i,:) =Ypos(1,:);
    
    %plotting here
    for j = 1:1:displaySize
        plot(prevX(:,j),prevY(:,j),'color',[colour(1,j) 0 j/displaySize])
        xlim([0 workX])
        ylim([0 workY])
        legend(['Temperature:' num2str(averageTemp)])
        drawnow
        hold on
    end
    
    
end