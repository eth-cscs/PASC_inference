%clear all

%% generate some funny data
%X_exact = ones(100,150);
X_exact = ones(20,30,10);

disp('- preparing movie')
tic
for t=1:size(X_exact,3)
    c = [size(X_exact,1),size(X_exact,2)].*[0.3+t*0.05,0.5+t*0.01];
%    c = [size(X_exact,1),size(X_exact,2)].*[0.3,0.5];
    
    r = size(X_exact,1)*0.3;
    for i=1:size(X_exact,1)
        for j=1:size(X_exact,2)
            if (i-c(1))^2+(j-c(2))^2 <= r^2
                X_exact(i,j,t) = 0;
            end
        end
    end
end
X = X_exact + 5*(-0.5 + rand(size(X_exact))); % add noise
X = min(max(X,0),1); % cut data
mytime = toc;
disp(['  finished in ' num2str(mytime) 's'])

%% solve K-means problem
% set number of clusters and penalty
T = size(X,3);
R = size(X,1)*size(X,2);
K = 2;
epssqr = 10^0;
DG = get_graph_matrix_grid(size(X,1),size(X,2),size(X,3),0.5);
theta_given = [0,1];
[ X_rec, theta, gamma, it, Lit ] = kmeans_H1_graph(X,K,DG,epssqr,theta_given);

% plot data
figure
for t=1:T
    ax = subplot(3,T,t);
    hold on
    title(['t = ' num2str(t)])
    myimage_show(:,:,1) = X_exact(:,:,t);
    myimage_show(:,:,2) = X_exact(:,:,t);
    myimage_show(:,:,3) = X_exact(:,:,t);
    image(min(max(myimage_show,0),1));
    axis(ax,'image');
    set(ax,'XTickLabel','','YTickLabel','')
    hold off
    
    ax = subplot(3,T,t+T);
    hold on
    myimage_show(:,:,1) = X(:,:,t);
    myimage_show(:,:,2) = X(:,:,t);
    myimage_show(:,:,3) = X(:,:,t);
    image(min(max(myimage_show,0),1));
    axis(ax,'image');
    set(ax,'XTickLabel','','YTickLabel','')
    hold off
    
    ax = subplot(3,T,t+2*T);
    hold on
    myimage_show(:,:,1) = X_rec(:,:,t);
    myimage_show(:,:,2) = X_rec(:,:,t);
    myimage_show(:,:,3) = X_rec(:,:,t);
    image(min(max(myimage_show,0),1));
    axis(ax,'image');
    set(ax,'XTickLabel','','YTickLabel','')
    hold off
end

% plot gamma
R = size(X,1)*size(X,2);
figure
title('Model indicator functions')
for t = 1:T
    for k = 1:K
        ax = subplot(K,T,(k-1)*T+t);
        hold on
        if k==1
            title(['t = ' num2str(t)])
        end
        Gammak = reshape(gamma((k-1)*R*T+t : T : k*R*T),size(X,1),size(X,2));
        myimage_show(:,:,1) = Gammak;
        myimage_show(:,:,2) = Gammak;
        myimage_show(:,:,3) = Gammak;
        image(min(max(myimage_show,0),1));
        axis(ax,'image');
        set(ax,'XTickLabel','','YTickLabel','')
        hold off
    end
end
     
% plot decrease of object function value (should be monotonically decreasing)
if isempty(theta_given)
    figure
    hold on
    title('Objective function during iterations')
    plot(1:length(Lit),Lit,'r')
    xlabel('$it$','Interpreter','latex')
    ylabel('$L_{it}$','Interpreter','latex')
    hold off
end
