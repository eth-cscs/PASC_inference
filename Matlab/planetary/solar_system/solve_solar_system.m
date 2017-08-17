function [t,y,bodies,dat]=solve_solar_system(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to solve for the motion of the planets in our solar system.    %
%                                                                         %
% Written by Dr Paul Connolly 2015 (copyright 2015).                      %
% Code provided as-is without any warranty.                               %
%                                                                         %
% Run by typing:                                                          %
% [t,y,bodies,dat]=solve_solar_system(1:10,trun);                         %
% where:                                                                  %
% 1:10 indicates to include interactions with all 10 bodies               %
% trun is the time in years you want the model to run for.                %
%                                                                         %
% Or by typing:                                                           %
% [t,y,bodies,dat]=solve_solar_system(1,trun);                            %
% where:                                                                  %
% 1 indicates to include interactions with the sun only.                  %
% trun is the time in years you want the model to run for.                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% some variables to use (don't change these)+++++++++++++++++++++++++++++++
global n_bodies;
global G;
global m;
global Gm;
global interaction;
global interactions;
global time1;
global sign1;
global starti;
%--------------------------------------------------------------------------


% Stolen by star: n_bodies=11; mstar=5e30; zstar=1e9*1e3;
% Binary star system: n_bodies=11; mstar=m(1)./2; zstar=1e9.*1e3; uxstar=sqrt(F*zstar/(m));


% SET THE NUMBER OF BODIES ++++++++++++++++++++++++++++++++++++++++++++++++
% The order is: Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus,
% Neptune and Pluto, a total of 10 bodies
n_bodies=10; % This can be set to less than 10 
if(nargin>0)
    n_bodies=varargin{1};
end
%--------------------------------------------------------------------------


% SET THE BODIES THAT A SINGLE BODY INTERACTS WITH (INCLUDING ITSELF!)+++++
% AND A FLAG FOR THE BODIES TO CONSIDER +++++++++++++++++++++++++++++++++++
inds1=ones(n_bodies,1); %consider all bodies in the solar system
% inds1=[1 0 0 1 0 1 0 0 0 0 0]; % consider the sun, earth and jupiter
% inds1=[1 0 0 0 0 0 0 0 0 0 1]; % sun and other star (n_bodies=11)

% interactions should index the bodies being simulated
interactions=1:n_bodies; % this means all bodies interact with each other. 
% interactions=1; % this means all bodies interact with body no 1, which is the sun. 
if(nargin>1)
    interactions=varargin{2};
end
if(nargin>2)
    inds1=varargin{3};
end
%--------------------------------------------------------------------------


% SET THE TIME ARRAY TO GET SOLUTION ON IN YEARS ++++++++++++++++++++++++++
% Time to get solution on in years (default - set again based on input)++++
% tsol=[0:-0.01:-200]; 
tsol=[0:0.01:1000]; 
if(nargin>4)
    tsol=[0:0.01:varargin{5}];
end
%--------------------------------------------------------------------------

starti=1;
if nargin>3
    sun_flag=varargin{4};
    if sun_flag == 0
        starti=2;
    end
end


inds=find(inds1==1);


G=6.67e-11;
m=[1.989e30 5e30]; % the mass of the sun and other star (n_bodies=11)
zstar=1e12; % distance of other star away from sun (n_bodies=11)
%uxstar=sqrt(G.*m(1)./zstar);
uxstar=0.;


% INITIAL DATA FOR THE SOLAR SYSTEM (TAKEN FROM JPL EPHEMERIS) ++++++++++++
% the product of G and m for the bodies in the solar system
Gm=[G.*m(1)./1e9 22032.09 324858.63 398600.440 ...
    42828.3 126686511 37931207.8 ...
    5793966 6835107 872.4 G.*m(2)./1e9].*1e9;

% The positions (x,y,z) and the velocities (vx,vy,vz) of all the planets
x=[0 1.563021412664830E+07 -9.030189258080004E+07 -1.018974476358996E+08 ...
    -2.443763125844157E+08 -2.35165468275322006E+08 -1.011712827283427E+09 ...
    2.934840841770302E+09 4.055112581124043E+09 9.514009594170194E+08 0].*1e3;
y=[0 4.327888220902108E+07 5.802615456116644E+07 1.065689158175689E+08 ...
    4.473211564076996E+07 7.421837640432589E+08 -1.077496255617324E+09 ...
    6.048399137411513E+08 -1.914578873112663E+09 -4.776029500570151E+09 0].*1e3;
z=[0 2.102123103174893E+06 6.006513603716755E+06 -3.381951053601424E+03 ...
    6.935657388967808E+06 2.179850895804323E+06 5.901251900068215E+07 ...
    -3.576451387567792E+07 -5.400973716179796E+07 2.358627841705075E+08 zstar./1e3].*1e3;

ux=[0 -5.557001175482630E+01 -1.907374632532257E+01 -2.201749257051057E+01 ...
    -3.456935754608896E+00 -1.262559929908801E+01 6.507898648442419E+00 ...
    -1.433852081777671E+00 2.275119229131818E+00 5.431808363374300E+00 uxstar./1e3].*1e3;
uy=[0 1.840863017229157E+01 -2.963461693326599E+01 -2.071074857788741E+01 ...
    -2.176307370133160E+01 -3.332552395475581E+00 -6.640809674126991E+00 ...
    6.347897341634990E+00 4.942356914027413E+00 -2.387056445508962E-02 0].*1e3;
uz=[0 6.602621285552567E+00 6.946391255404438E-01 1.575245213712245E-03 ...
    -3.711433859326417E-01 2.962741332356101E-01 -1.434198106014633E-01 ...
    4.228261484335974E-02 -1.548950389954096E-01 -1.551877289694926E+00 0].*1e3;
%--------------------------------------------------------------------------



% descriptor of the body being modelled++++++++++++++++++++++++++++++++++++
bodies={'Sun','Mercury','Venus','Earth','Mars','Jupiter','Saturn',...
    'Uranus','Neptune','Pluto','Black hole'};


% Mean distance from the sun+++++++++++++++++++++++++++++++++++++++++++++++
meandist=[7e8 5.79e10 1.082e11 1.496e11 2.279e11 7.783e11 1.426e12 ...
    2.871e12 4.497e12 5.914e12 sqrt(x(end).^2+y(end).^2+z(end).^2)];

% This indexes the bodies by what we want to consider++++++++++++++++++++++
Gm=Gm(inds);
x=x(inds);
y=y(inds);
z=z(inds);
ux=ux(inds);
uy=uy(inds);
uz=uz(inds);
bodies=bodies(inds);
meandist=meandist(inds);

n_bodies=length(inds);
%--------------------------------------------------------------------------


% set x,y,z,vx,vy,vz for each consecutive planet and the tolerances for the
% solution
YINIT=[];
AbsTol=[];
for i=1:n_bodies
    YINIT=[YINIT x(i) y(i) z(i) ...
        ux(i) uy(i) uz(i)];
    if(i==1)
        AbsTol=[AbsTol [1./1e3
                1./1e3
                1./1e3
                1./1e6
                1./1e6
                1./1e6]'];
    else
        AbsTol=[AbsTol [sqrt(x(i).^2+y(i).^2+z(i).^2)./meandist(i)./1e6
                sqrt(x(i).^2+y(i).^2+z(i).^2)./meandist(i)./1e6
                sqrt(x(i).^2+y(i).^2+z(i).^2)./meandist(i)./1e6
                sqrt(ux(i).^2+uy(i).^2+uz(i).^2)./1e7 
                sqrt(ux(i).^2+uy(i).^2+uz(i).^2)./1e7
                sqrt(ux(i).^2+uy(i).^2+uz(i).^2)./1e7]'];
    end
end
if(n_bodies==11) % for the black-hole, which has initial velocity zero (if simulated)
    AbsTol(end-5:end)=AbsTol(1:6)./100;
end
%--------------------------------------------------------------------------

% Pass options to the solver routine+++++++++++++++++++++++++++++++++++++++
sign1=sign(tsol(end)-tsol(1));
options=odeset('RelTol',1e-6,'AbsTol',AbsTol,'OutputFcn',@myplotfun); 
time1=0.;
% -------------------------------------------------------------------------


% Solve this problem using ODE113 +++++++++++++++++++++++++++++++++++++++++
[t,y]=ode113(@solar01,tsol.*365.25.*24.*3600,YINIT,options);
% -------------------------------------------------------------------------


% now parse y into something more human readable:++++++++++++++++++++++++++
for i=1:n_bodies
    eval(['dat.',bodies{i},'.xyz=y(:,[1:3]+(i-1).*6);']);
    eval(['dat.',bodies{i},'.uvw=y(:,[4:6]+(i-1).*6);']);
end
%--------------------------------------------------------------------------

% Function describing the time derivatives to be integrated++++++++++++++++
function dy=solar01(t,y)

global n_bodies;
global G;
global m;
global Gm;
global interaction;
global interactions;
global starti;

dy=zeros(6.*(n_bodies),1);
% x, y, z, vx, vy, vz
for i=starti:n_bodies
    ind=interactions; % bodies we are considering interactions between
    i2=find(ind==i);  % a body can't interact with itself
    ind(i2)=[];

    % square of distance between this planet and the other objects
    R2=(y((i-1).*6+1)-y((ind-1).*6+1)).^2+...
        (y((i-1).*6+2)-y((ind-1).*6+2)).^2+...
        (y((i-1).*6+3)-y((ind-1).*6+3)).^2;
    % inverse square law between body i and the rest of them
    dy((i-1).*6+4)=dy((i-1).*6+4) ...
        -sum(Gm(ind)'./R2.*(y((i-1).*6+1)-y((ind-1).*6+1))./sqrt(R2));% dvx/dt
    dy((i-1).*6+5)=dy((i-1).*6+5) ...
        -sum(Gm(ind)'./R2.*(y((i-1).*6+2)-y((ind-1).*6+2))./sqrt(R2));% dvy/dt
    dy((i-1).*6+6)=dy((i-1).*6+6) ...
        -sum(Gm(ind)'./R2.*(y((i-1).*6+3)-y((ind-1).*6+3))./sqrt(R2));% dvz/dt

    % rate of change of position, because we are solving a second order ODE
    % - always the same
    dy((i-1).*6+1)=y((i-1).*6+4);
    dy((i-1).*6+2)=y((i-1).*6+5);
    dy((i-1).*6+3)=y((i-1).*6+6);
end
% -------------------------------------------------------------------------


% Function to output time++++++++++++++++++++++++++++++++++++++++++++++++++
function status = myplotfun(t,y,flag)
global time1;
global sign1;
if(sign1==1) % check if the solution is going forward or backward-in-time
    if(t./86400./365.25>time1+1) % if last output was 1 year ago
        time1=time1+1;
        disp(['Time: ',num2str(t./86400./365.25),' yrs']);
    end
else
    if(t./86400./365.25<time1-1) % if last output was 1 year in future 
        time1=time1-1;
        disp(['Time: ',num2str(t./86400./365.25),' yrs']);
    end
end    
status=0.;
%--------------------------------------------------------------------------

