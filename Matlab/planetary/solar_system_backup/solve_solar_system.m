function [t,y,bodies]=solve_solar_system(planets, interactions_in, run_time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to solve for the motion of the planets in our solar system.    %
% 
% originally written by Dr Paul Connolly 2015 (copyright 2015).           %
% Code provided as-is without any warranty.                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% some variables to use (don't change these)+++++++++++++++++++++++++++++++
global n_bodies;
global G;
global m;
global Gm;
global interactions;
global time1;
global sign1;
global starti;
%--------------------------------------------------------------------------

interactions = interactions_in;

% SET THE NUMBER OF BODIES ++++++++++++++++++++++++++++++++++++++++++++++++
n_bodies=planets.nmb; 

% SET THE BODIES THAT A SINGLE BODY INTERACTS WITH (INCLUDING ITSELF!)+++++
% AND A FLAG FOR THE BODIES TO CONSIDER +++++++++++++++++++++++++++++++++++
inds1=ones(n_bodies,1); %consider all bodies in the solar system
%--------------------------------------------------------------------------


% SET THE TIME ARRAY TO GET SOLUTION ON IN YEARS ++++++++++++++++++++++++++
% Time to get solution on in years (default - set again based on input)++++
% tsol=[0:-0.01:-200]; 
tsol=[0:0.01:run_time]; 
%--------------------------------------------------------------------------

starti=1;

inds=find(inds1==1);


G=6.67e-11;
m=1.989e30; % the mass of the sun

% INITIAL DATA FOR THE SOLAR SYSTEM (TAKEN FROM JPL EPHEMERIS) ++++++++++++
% the product of G and m for the bodies in the solar system
Gm=[G.*m./1e9 22032.09 324858.63 398600.440 ...
    42828.3 126686511 37931207.8 ...
    5793966 6835107 872.4].*1e9;

% The positions (x,y,z) and the velocities (vx,vy,vz) of all the planets
x=planets.x;
y=planets.y;
z=planets.z;

ux=planets.ux;
uy=planets.uy;
uz=planets.uz;
%--------------------------------------------------------------------------


% descriptor of the body being modelled++++++++++++++++++++++++++++++++++++
bodies=planets.names;

% Mean distance from the sun+++++++++++++++++++++++++++++++++++++++++++++++
meandist=planets.meandist;

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

