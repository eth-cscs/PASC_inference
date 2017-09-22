% Script Solution for 3-body problem (earth, moon, spaceship)
% 2-D
% Spaceship thruster is fired for a certain period of time from a...
% geostationary orbit
% Figure eight solution

clear;
close all;

global me mm ms G xm thrust

% -----Constants------%
            me = 59742e24;  % mass of the earth [kg]
            mm = 7.35e22;   % mass of the moon [kg]
            ms = 1.5e5;     % mass of the spaceship [kg]
            G  = 6.673e-11; % Gravitational Constant [m^3 kg^-1 s^-2]
            xm = 384e6;     % distance from the earth to the moon [m]

% -----Parameters------%
thrustduration = 200;   % [s]
        tfinal = 5e4;   % [s]

            x0 =  1; % initial position in x
            y0 =  1; % initial position in y
            v0 = sqrt((G*me)/(6.4e6+2e5));    % initial velocity
           vx0 =  1; % initial velocity in x
           vy0 =  1; % initial velocity in y

             t =  1; % time
      t_thrust =  1; % time for initiate thrust?
      t_travel =  1; % time for spaceship to travel
         tsize =  1; % time size?
         tspan =  [0,tfinal] ; % time span, entire?
   
%----If Statement for Thrust----% 
        thrust = 1; % amount of thrust [lb]

            y0 = [x0,y0,vx0,vy0]; % not sure what all needs to be here

         [t,y] = ode45('space_ship_programODE', tspan, y0);

             x = y(:,1) ; % position in x after time t
             y = y(:,2) ; % position in y after time t
           vx1 = y(:,3) ; % velocity in x after time t
           vy1 = y(:,4) ; % velocity in y after time t

%----Plotting Velocity Analysis----%
figure(1);
subplot(2,1,1), plot(t,vx1);
subplot(2,1,2), plot(t,vy1);
xlabel('Time (sec.)');
ylabel('Velocity (m/sec.)');


%----Plotting Earth and Moon----% 

%Earth%
figure(2);
h=0; k=0; r=6.4e6; N=256;
t=(0:N)*2*pi/N;
plot( r*cos(t)+h, r*sin(t)+k);
axis('equal')
hold on

%Moon%
h=xm; k=0; r=1.75e6; N=256;
t=(0:N)*2*pi/N;
plot( r*cos(t)+h, r*sin(t)+k);
axis('equal')

%----Plotting Orbit----%
plot(x,y);
title('Spaceship Orbit');
ylabel('y(t)');
xlabel('x(t)');
axis('equal')
hold off