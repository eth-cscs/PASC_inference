function ydot = space_ship_programODE(t,var)

global me mm ms G xm thrust

x  = var(1);    % x position [m]
y  = var(2);    % y position [m]
vx = var(3);    % x velocity [m/s]
vy = var(4);    % y velocity [m/s]

% -----2-D problem => 4 ODE's------%

yd1 = vx;
yd2 = vy;
yd3 = 1/ms*(-G*me*ms/(x^2+y^2)*x/sqrt(x^2+y^2) + G*mm*ms/((xm-x)^2+y^2)*x/sqrt((xm-x)^2+y^2) + thrust*x/sqrt(x^2+y^2)); % unit vectors instead of angles
yd4 = 1/ms*(-G*me*ms/(x^2+y^2)*y/sqrt(x^2+y^2) - G*mm*ms/((xm-x)^2+y^2)*y/sqrt((xm-x)^2+y^2) + thrust*y/sqrt(x^2+y^2)); % unit vectors instead of angles

ydot = [yd1;yd2;yd3;yd4];

end