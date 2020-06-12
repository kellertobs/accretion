
% set parameters
N  = 10;   % number of bodies
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

for k = 1:K
    
    % plot model progress
    if ~mod(k-1,25)  % <- lower number for more frequent plots
        figure(1); clf;
        scatter(X(1,1),X(1,2),M(1).*100,(V(1,1).^2+V(1,2).^2).^0.5,'filled'); hold on; colorbar; axis equal;
        scatter(X(2:N,1),X(2:N,2),M(2:N).*2e5,(V(2:N,1).^2+V(2:N,2).^2).^0.5,'filled');
        axis([-20 20 -20 20]);
        drawnow;
    end
    
    % loop through all bodies
    for ni = 1:N
        for nj = 1:N
            R       = sum((X(ni,:)-X(nj,:)).^2).^0.5 + 1e-32;  % compute distance between ni, nj
            F       = - G*M(ni)*M(nj)/R^2 .* (X(ni,:)-X(nj,:))./R;  % computing force between ni, nj
            V(ni,:) = V(ni,:) + F/M(ni) * dt;  % update velocity of body ni
        end
    end
    
    X = X + V.*dt;  % update position of all bodies
    
end