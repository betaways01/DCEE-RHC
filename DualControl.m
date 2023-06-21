%%
% The agent uses a combination of a particle filter and
% a Receding Horizon Control (RHC) algorithm to make decisions about its movements.
%%
close all
clear
clc

addpath('Functions')

UAVVel=2;
sampleTime=10;
bLim=900;
% Simulated source parameters
% true source
s.Q = 20; % Release rate per t
% source coodinates
s.x = 5;
s.y = 8;
s.z= 4;
s.u = 10;
s.phi = 1;
s.ci = 0.1;%0.14;  % Also s.D for Pasquil model
s.cii = 0.6;%0.53; % Also s.t or tau for Pasquil Model
% s.ci = 0.14;
% s.cii = 0.53;
% s.duration = 0;
m.thresh = 5e-4;

% Create rectangular domain area
xmin = 0;
xmax = 10;
ymin = 0;
ymax = 10;
zmin = 0;
zmax = 5;
domain = [xmin xmax ymin ymax]; % Size of search area

% example data
stepsize = 0.1; %horisontal (x and y) spacing

x_coord = xmin : stepsize : xmax;
y_coord = ymin : stepsize : ymax;
z_coord = zmin : stepsize : zmax;

% Create 3D grid
[X,Y,Z] = meshgrid(x_coord,y_coord,z_coord);

ex.x_matrix = X;
ex.y_matrix = Y;
ex.z_matrix = Z;

% Initialisation and parameters
StartingPosition = [2 2,4]; % Starting position [x,y,z]
moveDist = 0.5; % How far to move


P_k = StartingPosition; % Current position
P_k_store = [];
P_k_store = [P_k_store; P_k];

pos.x_matrix = P_k(1);
pos.y_matrix = P_k(2);
pos.z_matrix = P_k(3);

% Plot example dispersion from true source
% conc = simpleGaussianPlume(s,m,ex);
conc = RadioactiveDispersionModel(s,ex);
figure(1)
hold off
height = 1;
concSurf=conc(:,:,height)/s.Q;
concSurf(concSurf<=m.thresh)=NaN;
pcolor(ex.x_matrix(:,:,height),ex.y_matrix(:,:,height),concSurf);
axis([xmin xmax ymin ymax]);
% slice(conc,28,17,1);
shading interp
% axis equal
% view(0,90)
xlab = xlabel('x');
ylab = ylabel('y');
colorbar
hold on
plot(s.x,s.y,'k.','MarkerSize',30)
plot(P_k(1),P_k(2),'rx')
plot(P_k_store(:,1),P_k_store(:,2),'r')
S = [];
colorbar

% initialise PF
N = 20000; %20000
PF_Memory=10;
resample = 0;
% Uniform prior for location
theta.x = xmin + (xmax-xmin) * rand(N,1);
theta.y = ymin + (ymax-ymin) * rand(N,1);
theta.z = ones(N,1)*s.z;
a = ones(N,1)*2;
b = ones(N,1)*5;
theta.Q = gamrnd(a,b);%200*rand(N,1);%
figure;histogram(theta.Q)
theta.u =s.u*ones(N,1);%2+6*rand(N,1);%0.75+0.5*rand(N,1);0 + randn(N,1)*0.5;%
theta.phi = s.phi*ones(N,1);%(10 + 30*rand(N,1)).*pi/180;
theta.ci = s.ci*ones(N,1);%0.12+0.1*rand(N,1);
theta.cii =s.cii*ones(N,1);%0.5+ 0.1*rand(N,1);
%Wp refers to particle weights
Wp = ones(N,1);
Wpnorm = Wp./sum(Wp);
Wp = Wpnorm;
  
figure(4)
hold on
% preprocess(s,theta,Wpnorm);

timestamp(1)=0;
D=[];
sampleHistory=P_k(1:2);

% Initialize arrays to store estimated source locations
estimated_x = [];
estimated_y = [];

% Main Loop
for i = 1:100
    % simulated data
    %Dsim = simpleGaussianPlume(s,m,pos);
    Dsim = RadioactiveDispersionModel(s,pos);
    % Dsim(Dsim<m.thresh)=0;
    ersize = size(Dsim);
    error = 0.05 * Dsim .* randn(ersize(1),ersize(2));
    Dsim = Dsim + error;
    Dsim(Dsim<m.thresh)=0;
    
    if rand<0.3
        Dsim = 0;
    end
    D(i)=Dsim;
    thetaPrev=theta;
    [theta, Wpnorm]=UpdatePFPlume(D, theta, Wpnorm, pos,P_k_store,m,N,PF_Memory,domain);

    % Call the RHC function to get the next position
    nextPosition = RHC(P_k);
    % Update the position
    P_k = nextPosition;
    RMSE_hist(i)=sqrt(mean(vecnorm([theta.x,theta.y]'-[s.x,s.y]')'.^2));
    
    figure(1)
    hold off
    concSurf=conc(:,:,height)/s.Q;
    concSurf(concSurf<=0.045)=NaN;
    pcolor(ex.x_matrix(:,:,1),ex.y_matrix(:,:,1),log10(concSurf))
    c1=min(min(concSurf));
    c2=max(max(concSurf));
    clim([-1.3 log10(c2)]);
    c=colorbar;
    c.Label.String='Concentration log(kg/m^3)';
    c.Limits=[-1.3,0];
    shading interp
    grid on
    hold on
    S(i) = 5+ceil(D(i)*5e4);
    scatter3(theta.x,theta.y,theta.z,3,'g','filled')
    plot3(pos.x_matrix,pos.y_matrix,pos.z_matrix,'ro','MarkerFaceColor','r','MarkerSize',5)
    plot3(P_k_store(:,1),P_k_store(:,2),P_k_store(:,3),'r--','LineWidth',2)
    plot3(s.x,s.y,s.z,'k.','markersize',20)
    xlab = xlabel('x (m)');
    ylab = ylabel('y (m)');
    set(xlab,'FontSize',16);
    set(ylab,'FontSize',16);
    set(gca,'fontsize',16)

    view(0,90)
    
    axis([xmin xmax ymin ymax])
    drawnow
    
%     if bLim<=810
%         if bLim<=630
%             if bLim<=400
%                 if bLim<=0
%                     true
%                 end
%             end
%         end
%     end
        
    
    indx = resampleStratified(Wpnorm);
    t.x = theta.x(indx);
    t.y = theta.y(indx);
    
    Xneighbour = zeros(1,5);
    Yneighbour = zeros(1,5);
    Zneighbour = zeros(1,5);
    
    %All proposed sensor locations are right,up,down,left for 1 step, 2
    %step and 3 step
    ynew = [[moveDist,moveDist,0,-moveDist,-moveDist,-moveDist,0,moveDist]];% ...
    %2*[0,moveDist,-moveDist,0] ...
    %3*[0,moveDist,-moveDist,0]];
    xnew = [[0,moveDist,moveDist,moveDist,0,-moveDist,-moveDist,-moveDist]];% ...
    %2*[moveDist,0,0,-moveDist] ...
    %3*[moveDist,0,0,-moveDist]];
    znew = [0,0,0,0,0,0,0,0,0,0,0,0];
    tic
    
    M= 40;%16;%5
    MM = 1;
    
    % Entrotaxis reward
    indx = resampleStratified(Wpnorm,M);
    d.x = theta.x(indx);
    d.y = theta.y(indx);
    d.z = theta.z(indx);
    d.Q = theta.Q(indx);
    d.u = theta.u(indx);
    d.phi = theta.phi(indx);
    d.ci = theta.ci(indx);
    d.cii = theta.cii(indx);
    %         d.duration = 0;%theta.duration(indx);
    
    
    
    %% ---------Entropy Reduction prediction at potential locations
    for k = 1:8 %for all proposed locations
        Xneighbour(k) = pos.x_matrix+xnew(k);
        Yneighbour(k) = pos.y_matrix+ynew(k);
        Zneighbour(k) = pos.z_matrix+znew(k);
        if pos.x_matrix+xnew(k)<xmin || pos.x_matrix+xnew(k)>xmax || pos.y_matrix+ynew(k)<ymin || pos.y_matrix+ynew(k)>ymax || pos.z_matrix+znew(k)<zmin || pos.z_matrix+znew(k)>zmax% || pos.x_matrix+xnew(k)>x_distance
            var(k)=NaN;
            dist_theta(k)=NaN;
            theta_RMSE(k)=NaN;
            continue
            %                     elseif any(all([Xneighbour(k),Yneighbour(k)]'==sampleHistory'))
            %                         var(k)=NaN;
            %                         dist_theta(k)=NaN;
            %                         theta_RMSE(k)=NaN;
            %                         continue
        end
        
        
        npos.x_matrix = Xneighbour(k);
        npos.y_matrix = Yneighbour(k);
        npos.z_matrix = Zneighbour(k);
        
        %pC = simpleGaussianPlume(theta,m,npos);
        pC = RadioactiveDispersionModel(theta,npos);
        
        %         pC(pC<m.thresh)=0;
        entropy=0;
        theta_RMSE(k)=0;
        dist_theta(k)=0;
        
        %desC = simpleGaussianPlume(d,m,npos);
        desC = RadioactiveDispersionModel(d,npos);
        designd = repmat(desC,1,MM);
        designd = designd+(1*designd.*randn(M,MM));
        %desC(desC<m.thresh)=0;
        %missdet = rand(M,MM);
        designd(rand<0.3)=0;
        designd(designd<m.thresh)=0;
        for jj = 1:M
            for jjj = 1:MM
                
                dC = designd(jj,jjj);
                
                
                zWp= Likelihood_Like_Yee(pC, dC, Wpnorm,m);
                zWpnorm = zWp./sum(zWp);

                theta_mean_xy=N*[mean(theta.x.*zWpnorm),mean(theta.y.*zWpnorm)];
                theta_RMSE(k)=theta_RMSE(k)+sum(zWpnorm.*vecnorm([theta.x,theta.y]'-theta_mean_xy')'.^2)/(M*MM);              
                dist_theta(k)=dist_theta(k)+norm(theta_mean_xy-[Xneighbour(k),Yneighbour(k)])/(M*MM);

            end
        end
        % Entrotaxis
        %         p_d = p_d./sum(p_d);
        %         entropy = entropy - (sum(log2(p_d+(p_d==0))));
        
        
        
        var(k) = dist_theta(k)+theta_RMSE(k);
        %var(k) = theta_RMSE(k);
        %var(k) = dist_theta(k);
        
    end
    
    %%
    %         var
    % if D(end)>m.thresh
    %     var([2 6 10])=var([2 6 10])*0.8;
    % end
    [val,ind] = min(var);
    [val2,ind2] = min(theta_RMSE);
    [val3,ind3] = min(dist_theta);
    dualControlJ(i,:)=[dist_theta(ind),theta_RMSE(ind),ind3==ind,ind2==ind];
    
    
    pos.x_matrix = Xneighbour(ind);
    pos.y_matrix = Yneighbour(ind);
    pos.z_matrix = Zneighbour(ind);
    
    P_k = [pos.x_matrix pos.y_matrix pos.z_matrix];
    sampleHistory=[sampleHistory;P_k(1:2)];
    
    move_time=floor(norm(P_k-P_k_store(i,:))/UAVVel)+sampleTime;
    bLim=bLim-move_time
    if bLim<=0
        break
    end
    timestamp(i+1)=timestamp(i)+move_time;
    P_k_store = [P_k_store; P_k];
    
    
    
    Covar = cov(theta.x,theta.y);
    Spread = sqrt(Covar(1,1)+Covar(2,2))
    
    if Spread<0 % 3.5? 4? 5 default
        break
    end


    % Store estimated source locations at each time step
    estimated_x = [estimated_x; mean(theta.x)];
    estimated_y = [estimated_y; mean(theta.y)];
        
    
end


% 3D Trajectory of the Agent
figure(2)
plot3(P_k_store(:,1), P_k_store(:,2), P_k_store(:,3), 'r-');
hold on;
plot3(s.x, s.y, s.z, 'k.', 'MarkerSize', 20);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Trajectory of the Agent');
grid on;

% Concentration Slices
figure(3)
for height = 1:5
    subplot(1, 5, height);
    concSurf = conc(:,:,height)/s.Q;
    concSurf(concSurf <= m.thresh) = NaN;
    pcolor(ex.x_matrix(:,:,height), ex.y_matrix(:,:,height), concSurf);
    shading interp;
    title(['Height = ', num2str(height)]);
    colorbar;
end

% Histogram of Particle Weights
figure(4)
histogram(Wpnorm);
xlabel('Particle Weights');
ylabel('Frequency');
title('Histogram of Particle Weights');

% Heatmap of Sensor Readings
figure(5)
heatmap(D);
title('Heatmap of Sensor Readings');

% Agent's Estimated Source Location Over Time
figure(6)
plot(estimated_x, 'r');
hold on;
plot(estimated_y, 'b');
legend('Estimated X', 'Estimated Y');
xlabel('Time Step');
ylabel('Estimated Source Location');
title('Agent''s Estimated Source Location Over Time');

indx = resampleStratified(Wpnorm);

theta.x = theta.x(indx);
theta.y = theta.y(indx);
theta.z = theta.z(indx);
theta.Q = theta.Q(indx);
theta.u = theta.u(indx);
theta.phi = theta.phi(indx);
theta.ci = theta.ci(indx);
theta.cii = theta.cii(indx);
figure(4)
hold on
preprocess(s,theta,Wpnorm);


%% Define a new function for the RHC algorithm: 
% This function takes the current position of the agent as input 
% and return the next position of the agent.
% Inside this function, we implement the steps of the RHC algorithm.
function nextPosition = RHC(currentPosition)
    % Define the lookahead steps
    lookaheadSteps = 2;

    % Initialize the entropy values array
    entropyValues = zeros(1, 8);
    
    % Calculate the value and entropy for each candidate point
    for i = 1:8
        % Get the candidate point
        candidatePoint = getPointByIndex(currentPosition, i);
        
        % Calculate the value of the candidate point
        value = calculateValue(candidatePoint);
        
        % If the value is not satisfactory, continue to the next point
        if value < 0.4  % Define SOME_THRESHOLD as per your requirements
            continue;
        end
        
        % Predict the position after applying DCEE for lookaheadSteps
        predictedPosition = candidatePoint;
        for j = 1:lookaheadSteps
            predictedPosition = DCEE(predictedPosition);
        end
        
        % Calculate the entropy of the predicted position
        entropyValues(i) = calculateEntropy(predictedPosition);
    end
    
    % Check if all entropy values are zero, implying no satisfactory point was found
    if all(entropyValues == 0)
        error('RHC: No satisfactory point found.');
    end
    
    % Select the point with the highest entropy
    [~, index] = max(entropyValues);
    selectedPoint = getPointByIndex(currentPosition, index);
    
    % Move one step towards the selected point
    nextPosition = takeOneStep(currentPosition, selectedPoint);
end

%% This function assumes that the agent can move in 8 directions: 
% up, down, left, right, and the four diagonals. 
% The points are indexed from 1 to 8 in a clockwise manner 
% starting from the point above the current position.
function point = getPointByIndex(currentPosition, index)
    % Adapt this function to work with 3D points
    directions = [0 1 0; 1 1 0; 1 0 0; 1 -1 0; 0 -1 0; -1 -1 0; -1 0 0; -1 1 0];
    direction = directions(index, :);
    point = currentPosition + direction;
end

%%  the "value" of a point could be related to the likelihood 
% of the emitting source being at that point. 
% This likelihood could be calculated based on the information available
% to the agent, such as sensor readings.
function value = calculateValue(candidatePoint)
    % Calculate the value of the candidate point
    % This could be based on sensor readings, a model of the environment, etc.
    
    % For example, if you have a function called 'getSensorReading' that returns
    % the sensor reading at a given point, and the sensor reading is directly
    % proportional to the likelihood of the emitting source being at that point,
    % you could calculate the value as follows:
    
    sensorReading = getSensorReading(candidatePoint);
    
    % If the sensor reading is directly proportional to the value
    value = sensorReading;
    
    % If the sensor reading is not directly proportional to the value, you would
    % need to use a different method to calculate the value. This could involve
    % using a model of the environment, taking into account other factors, etc.
end

%% This function should take a position as input 
% (in this case, the predicted position) 
% and return the next position of the agent based on the DCEE algorithm.
function nextPosition = DCEE(predictedPosition)
    % Calculate the entropy of 8 points around the predicted position
    entropyValues = calculateEntropy(predictedPosition);
    
    % Select the point with the highest expected performance
    [~, index] = max(entropyValues);
    selectedPoint = getPointByIndex(predictedPosition, index);
    
    % Move one step towards the selected point
    nextPosition = takeOneStep(predictedPosition, selectedPoint);
end


%% Shanon Entropy
function entropy = calculateEntropy(predictedPosition)
    % Get the 8 points around the predicted position
    points = getPointsAround(predictedPosition);
    
    % Initialize entropy
    entropyValues = zeros(1, length(points));
    
    % Calculate the entropy for each point
    for i = 1:length(points)
        % Get the sensor readings for the point
        readings = getSensorReadings(points(i, :));
        
        % Calculate the entropy of the readings
        entropyValues(i) = -sum(readings .* log2(readings + (readings == 0)));
    end
    
    % Return the average entropy of the neighboring points
    entropy = mean(entropyValues);
end


%% the agent moves stepSize units in the direction of the selected point. 
% The direction is calculated as a unit vector pointing 
% from the current position to the selected point.
function newPosition = takeOneStep(currentPosition, selectedPoint)
    % Adapt this function to work with 3D points
    stepSize = 1;
    direction = selectedPoint - currentPosition;
    direction = direction / norm(direction);
    newPosition = currentPosition + stepSize * direction;
end

%%
% Implement the getSensorReading function
function reading = getSensorReading(point)
    % Example: Simulate a sensor reading based on a simple dispersion model
    % and the position of the sensor.
    
    % Parameters of the dispersion model (for example)
    sourceStrength = 10;
    dispersionCoefficient = 0.1;
    
    % Position of the source (for example) in 3D
    sourcePosition = [5, 5, 1]; % Example values
    
    % Calculate the distance from the point to the source in 3D
    distance = norm(point - sourcePosition);
    
    % Calculate the reading based on the dispersion model
    reading = sourceStrength / (2 * pi * dispersionCoefficient * distance^2);
end

% Implement the getSensorReadings function
function readings = getSensorReadings(points)
    % Example: Get sensor readings for multiple points
    
    % Initialize readings array
    readings = zeros(1, size(points, 1));
    
    % Loop through each point and get the sensor reading
    for i = 1:size(points, 1)
        readings(i) = getSensorReading(points(i, :));
    end
end

%%
function points = getPointsAround(centerPoint)
    % Define the 8 directions in which the agent can move in 3D space
    directions = [0 1 0; 1 1 0; 1 0 0; 1 -1 0; 0 -1 0; -1 -1 0; -1 0 0; -1 1 0];
    
    % Initialize the array to store the points
    points = zeros(8, 3);
    
    % Calculate the coordinates of the 8 points around the centerPoint
    for i = 1:8
        points(i, :) = centerPoint + directions(i, :);
    end
end
