% Script to run the model for with and without interactions between planets

% SECTION 0: DO NOT CHANGE THESE ++++++++++++++++++++++++++++++++++++++++++
INTERACT_WITH_SUN_ONLY=1;
INTERACT_WITH_ALL_BODIES=2;
INTERACT_WITH_SUN_AND_JUPITER=3;

MOVING_SUN=1;
FIXED_SUN=0;

ALL_BODIES=1;
JUST_EARTH_SUN_AND_JUPITER=2;
ALL_BODIES_EXCEPT_JUPITER=3;

FIRST=1;
SECOND=2;
%--------------------------------------------------------------------------


% SECTION 1: USER INPUTS - CHANGE THESE++++++++++++++++++++++++++++++++++++
which_bodies=ALL_BODIES;
which_interactions=INTERACT_WITH_SUN_ONLY;
which_star=MOVING_SUN;
which_output=FIRST;
run_time=10;
plot_figures=true;
plot_stride=10;
%--------------------------------------------------------------------------


% SECTION 2: PARSE USER INPUTS

% SELECT WHICH BODIES WE ARE CONSIDERING+++++++++++++++++++++++++++++++++++
switch which_bodies
    case ALL_BODIES
        bodies01=ones(1,10);
    case JUST_EARTH_SUN_AND_JUPITER
        bodies01=[1 0 0 1 0 1 0 0 0 0];
    case ALL_BODIES_EXCEPT_JUPITER
        bodies01=[1 1 1 1 1 0 1 1 1 1];        
    otherwise
        disp('error which_bodies');
        return;
end
%--------------------------------------------------------------------------


% SELECT WHICH BODIES INTERACT WITH WHICH++++++++++++++++++++++++++++++++++
switch which_interactions
    case INTERACT_WITH_SUN_ONLY
        inters01=1;
    case INTERACT_WITH_ALL_BODIES
        inters01=1:10;
    case INTERACT_WITH_SUN_AND_JUPITER
        inters01=[1 6];        
    otherwise
        disp('error which_interactions');
        return;
end
%--------------------------------------------------------------------------



% SECTION 3: RUN THE MODEL
switch which_output
    case FIRST
    [t,y,bodies,dat]=solve_solar_system(10,inters01,bodies01,...
        which_star,run_time); 
    case SECOND
    [t1,y1,bodies1,dat2]=solve_solar_system(10,inters01,bodies01,...
        which_star,run_time);
    otherwise
        disp('error which_output');
        return;        
end



