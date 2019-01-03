figure
h = animatedline;
axis([0 4*pi -1 1])

for x = linspace(0,4*pi,10000)
    y = sin(x);
    addpoints(h,x,y)
    pause(0.0001)
    drawnow limitrate
end