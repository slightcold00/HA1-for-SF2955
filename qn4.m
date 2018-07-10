load('stations.mat')
load('RSSI-measurements.mat')
N = 100000;
n = 501;
tau = zeros(2,n); % vector of estimates
v = (1.5^2);
trackW = zeros(N,n);



%f = @(x1,x2,y1,y2) pdist([x1,y1;x2,y2],'euclidean');
p = @(x,y) mvnpdf(y,[90-10*3*log10(cal(1,x,pos_vec));
                    90-10*3*log10(cal(2,x,pos_vec));
                    90-10*3*log10(cal(3,x,pos_vec));
                    90-10*3*log10(cal(4,x,pos_vec));
                    90-10*3*log10(cal(5,x,pos_vec));
                    90-10*3*log10(cal(6,x,pos_vec))]',diag([v,v,v,v,v,v],0)); % observation density, for weights

part = transpose(mvnrnd(zeros(6,1),diag([500,5,5,200,5,5],0),N));% initialization
w = p(part,Y(:,1)');
tau(1,1) = sum(part(1,:).*w')/sum(w);
tau(2,1) = sum(part(4,:).*w')/sum(w);
trackW(:,1) = w;

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
simulate_Z = simulate(mc,n);
for  k = 1:(n-1) % main loop
    zM = repmat(rZ*Z(:,simulate_Z(k)),1,N);
    wM = rW*transpose(mvnrnd([0;0],diag([0.25;0.25],0),N));
    xM = rX* part; 
    YmN = repmat(Y(:,k+1),1,N);
    part = xM + zM + wM;
    w =  p(part,Y(:,k+1)');
    %w = w/max(w);
    ind = randsample(N,N,true,w);
    part = part(:,ind);
    %for i = 1:N
    %    w(i) = w(i).*p(part(:,i),YmN(:,i));
    %end
    %w = w.*p(part,Y(:,k + 1)); % weighting
    %dummy = part.*w;
    tau(1,k+1) = sum(part(1,:).*w')/sum(w);
    tau(2,k+1) = sum(part(4,:).*w')/sum(w);
    %tau(k + 1) = sum(part.*w)/sum(w); % estimation
    %trackP(:,k+1) = part;
    %trackW(k+1,:) = w;
    trackW(:,k+1) = w;
end
%Plot points
figure
plot(tau(1,:),tau(2,:),'*'); hold on;
plot(pos_vec(1,:),pos_vec(2,:),'*','Color',[1 0 0]);
title('Simulated Trajectory')
xlabel('x1')
ylabel('x2')

% Calculating the efficency sampling
track_sample = zeros(1,n);
sample_size = 0;
dummy = Inf;
for i = 1:n
    CV2 = (1/N)*sum((N*(trackW(:,i)./sum(trackW(:,i)))-1).^2);
    track_sample(i) = N/(1+CV2);
    if track_sample(i)<dummy
        dummy = track_sample(i);
        sample_size = i;
       
    end
end
%histogram when n= 50
figure;
histogram(log(trackW(:,100)),20);
%histogram when n= 100
figure;
histogram(log(trackW(:,150)),20);
%histogram when n= 150
figure;
histogram(log(trackW(:,200)),20);
%efficient sample size
figure;
plot(1:n,track_sample)
%Plot histogram (Assuming that n=1)
%tt = track_Z(5,:);
%figure
%histogram(trackW(3,:))