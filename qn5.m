load('stations.mat')
load('RSSI-measurements-unknown-sigma.mat')
N = 100000;
n = 501;
tau = zeros(2,n); % vector of estimates
v = 0.5;
max_cN = log(0); %produces negative infinity
max_v = 0;
track_vm = zeros(1,13);
vm = 1;
trackW = zeros(N,n);
while v<=3
    vT = v^2;    
%f = @(x1,x2,y1,y2) pdist([x1,y1;x2,y2],'euclidean');
    p = @(x,y) mvnpdf(y,[90-10*3*log10(cal(1,x,pos_vec));
                    90-10*3*log10(cal(2,x,pos_vec));
                    90-10*3*log10(cal(3,x,pos_vec));
                    90-10*3*log10(cal(4,x,pos_vec));
                    90-10*3*log10(cal(5,x,pos_vec));
                    90-10*3*log10(cal(6,x,pos_vec))]',diag([vT,vT,vT,vT,vT,vT],0)); % observation density, for weights

    part = transpose(mvnrnd(zeros(6,1),diag(sqrt([500,5,5,200,5,5]),0),N));% initialization
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
    %Calculate likelihood
    cN = sum(log(w));
    %cN = (1/n)*(-(n+1)*log(N)+sum(log(sum(w))));
    for  k = 1:(n-1) % main loop
        zM = repmat(rZ*Z(:,simulate_Z(k)),1,N);
        wM = rW*transpose(mvnrnd([0;0],diag([0.25;0.25],0),N));
        xM = rX* part; 
        YmN = repmat(Y(:,k+1),1,N);
        part = xM + zM + wM;
        w =  p(part,Y(:,k+1)');
        ind = randsample(N,N,true,w);
        part = part(:,ind);
        tau(1,k+1) = sum(part(1,:).*w')/sum(w);
        tau(2,k+1) = sum(part(4,:).*w')/sum(w);
        trackW(:,k+1) = w;
        cN = cN + sum(log(w));
    end
    track_vm(1,vm) = cN;
    %cN = cN / (N^(n+1));
    vm=vm+1;
    if(cN > max_cN)
        max_cN = cN;
        max_v = v;
    end
    v = v + 0.2;
    cN = 0;
end
plot(1:13,track_vm);
%Plot points
%figure
%plot(tau(1,:),tau(2,:),'*'); 
%hold on;
%plot(pos_vec(1,:),pos_vec(2,:),'*','Color',[1 0 0]);

%Plot histogram (Assuming that n=1)
%tt = track_Z(5,:);
%figure
%histogram(trackZ(3,:))