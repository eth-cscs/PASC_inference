clear all
close all

inters=1; % 1=INTERACT_WITH_SUN_ONLY;[1 ... 1]=INTERACT_WITH_BODIES;
run_time=0.1;
plot_figures=true;
plot_stride=1;
t_step=0.0001;

% prepare planets
planets.nmb = 5;

planets.x=[0 1.563021412664830E+07 2.563021412664830E+07 -2.763021412664830E+07  2.763021412664830E+07].*1e3;
planets.y=[0 2.327888220902108E+07 3.327888220902108E+07 -3.327888220902108E+07 -3.327888220902108E+07].*1e3;
planets.z=[0 0 0 0 0].*1e3;

planets.ux=[0 -5.557001175482630E+01 -3.557001175482630E+01 3.857001175482630E+01 4.257001175482630E+01].*1e3;
planets.uy=[0 4.840863017229157E+01 4.840863017229157E+01 -4.840863017229157E+01 4.840863017229157E+01].*1e3;
planets.uz=[0 0 0 0 0].*1e3;

planets.names = {'Star','Planet1','Planet2','Planet3','Planet4'};

planets.meandist = [7e8 5.79e10 5.79e10 5.79e10 5.79e10];

% RUN THE MODEL
[t,y,bodies]=solve_solar_system(planets,inters,run_time,t_step); 

animate_solar_system2
