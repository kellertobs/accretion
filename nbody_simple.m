
% set parameters
N  = 2e3;   % number of bodies
K  = 1e6;   % number of time steps
dt = 0.01;  % size of time step
G  = 1;     % gravitational constant

% initialise body mass and position
M  = [2;rand(N-1,1).*0.0001];
X  = [0,0;randn(N-1,2).*5];

% calculate initial orbital velocity
r  = sum((X-X(1,:)).^2,2).^0.5 + 1e-32;
R  = (X-X(1,:))./r;
V  = sqrt(G.*M(1)./r) .* [R(:,2),-R(:,1)];


%Josh's comment
for k = 1:K
    
    % plot model progress
    if ~mod(k-1,25)  % <- lower number for more frequent plots
        figure(1); clf;
        plot(X(1,1),X(1,2),'ko','MarkerFaceColor','k','MarkerSize',M(1).^(1/3).*20); hold on; box on; axis equal;
        scatter(X(2:N,1),X(2:N,2),M(2:N).^(1/3).*2e2,sum(V(2:N,:).^2,2).^0.5,'filled'); caxis([0.25,1]); colorbar;
        axis([X(1,1)-20,X(1,1)+20,X(1,2)-20,X(1,2)+20]);
        title('Accretionary Disk');
        drawnow;
        figure(2); clf;
        histogram(log10(M),min(50,max(5,round(N/10))));
        title('Size Distribution');
    end
    
    % calculate new orbital velocity
    for nj = 1:N
        D  = sum((X-X(nj,:)).^2,2).^0.5 + 1e-16;
        Fj = - G.*(M.*M(nj))./D.^2 .* (X-X(nj,:))./D;
        V  = V + Fj./M .* dt;
    end
    
    X = X + V.*dt;  % update position of all bodies
    
    % detect collisions and merge collided bodies
    nj = 1;
    while nj < N
        D    = sum((X-X(nj,:)).^2,2).^0.5 + 1e-16;
        ind  = find(D < (M+M(nj)).^(1/3)./10,2);
        if length(ind)>1
            ind = ind(end);
            V(nj,:)   = (M(ind)*V(ind,:) + M(nj)*V(nj,:))/(M(ind)+M(nj));
            M(nj)     =  M(ind) + M(nj);
            X(ind,:)  =  [];
            V(ind,:)  =  [];
            M(ind)    =  [];
            N         =  N-1;
        end
        nj = nj+1;
    end
    
end