M = 10000;
T = 1; %final time in years in trajectory
dt = T/M; %time step
D = [500; 5; 5; 200; 5; 5];
X0 = transpose(mvnrnd(zeros(6,1),diag(sqrt([500,5,5,200,5,5]),0)));
P = 1/20*[16 1 1 1 1;1 16 1 1 1;1 1 16 1 1;1 1 1 16 1;1 1 1 1 16];
deltaT = 0.5;
alpha = 0.6;
% Matrices needed for Equation
rXC = [1 deltaT (deltaT^2)/2;0 1 deltaT;0 0 alpha];
rX = [rXC zeros(3,3); zeros(3,3) rXC];
rZC = [(deltaT^2)/2;deltaT;0];
rZ = [rZC zeros(3,1); zeros(3,1) rZC];
rWC = [(deltaT^2)/2;deltaT;1];
rW = [rWC zeros(3,1); zeros(3,1) rWC];
% Simulate Z
Z = [[0;0] [3.5;0] [0;3.5] [0;-3.5] [-3.5;0]];
mc = dtmc(P);
simulate_Z = simulate(mc,M);
% Create Initial Matrix
X = zeros(6,M+1);
X(:,1) = X0;
% Create simulation
for i=1:(M)
    zM = rZ*Z(:,simulate_Z(i));
    wM = rW*[normrnd(0,0.5); normrnd(0,0.5)];
    xM = rX* X(:,i); 
    X(:,i+1) = xM + zM + wM;
end
t = 0:dt:T;
figure
l= animatedline;
minLimit = min(min(X(1,:)),min(X(4,:)));
maxLimit = max(max(X(1,:)),max(X(4,:)));
axis([minLimit maxLimit minLimit maxLimit])
x1 = X(1,:);
x2 = X(4,:);
l2 = animatedline('Color','r');
l2.Marker = '*';
for i=1:length(x1)
    addpoints(l,x1(i),x2(i))
    %addpoints(l2,x1(i),x2(i))
    drawnow
    pause(0.01)
    
end
%g = hgtransform;
%patch('XData',X(:,1),'YData',X(:,2),'FaceColor','yellow','Parent',g)
%minLimit = min(min(X(:,1)),min(X(:,2)));
%maxLimit = max(max(X(:,1)),max(X(:,2)));
%axis equal
%xlim([minLimit maxLimit])
%ylim([minLimit maxLimit])
%for t=1:length(t)-1
 % dummy = [0];
 % pt1 = X(t,:);
 % pt1_D = [pt1 dummy];
  %pt2 = X(t+1,:);
 % pt2_D = [pt2 dummy];
 % g.Matrix = makehgtform('translate',pt1_D + t*(pt1_D-pt2_D));
 % drawnow
%end
%plot3(t,X(:,1),X(:,2));