%// Sample x and y values assumed for demo.
x = 1:1000;
y = x.^2;

%// Plot starts here
figure,hold on

%// Set x and y limits of the plot
xlim([min(x(:)) max(x(:))])
ylim([min(y(:)) max(y(:))])

%// Plot point by point
for k = 1:numel(x)
    plot(x(k),y(k),'-') %// Choose your own marker here

    %// MATLAB pauses for 0.001 sec before moving on to execue the next 
    %%// instruction and thus creating animation effect
    pause(0.001);     
end