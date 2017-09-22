mybackground = 0.51*ones(1,3);
myplanets = 0.5*ones(1,3);

%mybackground = 1.0*ones(1,3);
%myplanets = 0.0*ones(1,3);

% create animation of solar system
figure('renderer','openGL','position',[0    0   256   128],'Color',mybackground);
hold all
axis off
axis equal
axis([-5e10 5e10 -5e10 5e10])
view(0,90)
axis manual

if plot_figures==true
    mkdir planetary2/
end

% Plot the paths of the planets
if false
    for j=1:planets.nmb-1
        h2(j)=plot3(y(:,1+j*6)-y(:,1),y(:,2+j*6)-y(:,2),y(:,3+j*6)-y(:,3),'linewidth',0.5);
    end
end

xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
     
% Scale the sizes of the planets
sz_sun=70;
sz = 0.8*[40 40 35 38 42 44 30];
%sz=log([4850 12000 12000 6790 142980 120536 51118 49500 2368]./0.00001);
k=1;

% Loop over all frames
for i=1:plot_stride:length(t)
    hold on

    disp([num2str(i) '/' num2str(length(t))])
    
    %plot sun
    h_sun =plot3(0,0,0,...
           '.','color',myplanets,'markersize',sz_sun,'MarkerEdgeColor', myplanets,'MarkerFaceColor', myplanets);
    
    % plot planets
    for j=1:planets.nmb-1
        h(j)=plot3(y(i,1+j*6)-y(i,1),y(i,2+j*6)-y(i,2),y(i,3+j*6)-y(i,3),...
            '.','color',myplanets,'markersize',sz(j),'MarkerEdgeColor', myplanets,'MarkerFaceColor', myplanets);
    end
    pause(0.01);
    hold off
    if plot_figures==true
%        eval(['print -dpng pics/frame',num2str(k,'%03d'),'.png']);
%        fig = gcf;
%        set(gcf,'PaperPositionMode','auto')
%        set(gcf,'PaperUnits','normalized')
%        set(gcf,'PaperPosition',[0 0 400/1275 250/1650])
%        eval(['print -djpeg90 pics/frame',num2str(k,'%03d'),'.jpg']);

        FileName = ['planetary2/frame' num2str(k,'%04d')];

        img = getframe(gcf);
%        imwrite(img.cdata, [FileName, '.jpg'], 'Quality', 90);    
        imwrite(img.cdata, [FileName, '.bmp']);
    end
    k=k+1;
    delete(h);
    delete(h_sun);

end
