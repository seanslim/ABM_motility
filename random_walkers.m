%% TO DO

% x. track the total number of "stuck" vs number of "swimmers"
% 8. make some stuck forever
% x. histogram the number of times a cell is stuck
% 3. incorporate "energy taxis" where cells tumble more where there is a
% high pmf (tumbling rate is higher where glucose-g is high)
% 4. calculate steadistate of stuck vs unstuck


%% Setup
% units in micrometers
% speeds in micrometers / sec

close all
clearvars A A1 A2 c C d h1 h2 h3 Hist_stuck N0 n_stuck n_swim Nx Ny rms T V Vrec

tmax = 240;
tstep = 3;

iter = 1;
D = zeros(iter,3);

% size of area

x = 2000;
y = 2000;

for k = 1:iter
% number of cells

n = 200; %can go up to 10e5 cells without much problem
N0 = rand(n,2)*10+x/2;
Nx = N0(:,1);
Ny = N0(:,2);

%% Velocity distribution

v = [0 0 60 60 60 60 60 60 60 60]; % [0 0 45 45 60 60 60 60 75 75];

for i = 1:n
    V(i,1) = v(randi(length(v)));
end

Vrec = ones(n,tmax/tstep+2);
Vrec(:,1) = V(:,1);

%% Run time distribution and rules

rt = 3;
T = randi(rt)*ones(n,1); % needs to be altered by local conc of glucose, glucose needs to be used up by cells

%% Tumble distribution

A = rand(n,3)*360;
athresh = 150;

% [0 30 60 90 120 150 180];

%% Glucose field and consumption distribution

cbins = 20;
g = ones(cbins,cbins);

cons = 0.001;

%% Simulation

% figure() % rms evolution
rms = zeros(tmax/tstep+1,1);
% h2 = plot(0:tstep:tmax,rms);
% axis([0,tmax,0,1000])
% 
% h2.YDataSource = 'rms';
% 
figure() % scatter and color plot
c = [rand(n,1) rand(n,1) rand(n,1)];
nbins = 20;

X_h3 = x/(nbins*2):x/nbins:x-x/(nbins*2); %A2{1}; put below
Y_h3 = y/(nbins*2):y/nbins:y-y/(nbins*2); %A2{2};

[A1,A2] = hist3([Nx,Ny],'Edges',{X_h3 Y_h3});

% [C, h3] = contourf(A2{1},A2{2},A1);

h3 = pcolor(X_h3,Y_h3,A1);
hold on
h1 = scatter(N0(:,1),N0(:,2),25,c,'filled');
axis([0 x 0 y]);
hold off

h1.XDataSource = 'Nx';
h1.YDataSource = 'Ny';
% h3.XData = 'X_h3'; % for evolving window
% h3.YData = 'Y_h3';
% h3.CData = 'A1';

for t=0:tstep:tmax
    
%% Construct state vectors (stuck, moving, tumbling)
 
        % countdown
        % velocity
        % angle
    T= T-tstep;
    
    for i=1:n
        
        % A(:,1) = A(:,2);
        
        if T(i)<= 0 % if the run is over
            
            A(i,1) = A(i,2); % sets the history
            if V(i,1) ~= 0
                A(i,3) = A(i,2); % if its in the stuck state, it will maintain its A(i,3) value
            end
            A(i,2) = randi(360); % conditional change in direction
            
            % agar tunnel condition
            if V(i,1) ~= 0
            V(i,1) = v(randi(length(v))); % random velocity

            else % must have a certain angle to escape
            if abs(A(i,2)-A(i,3)) >= athresh && abs(A(i,2)-A(i,3)) <= athresh+(180-athresh)*2
                V(i,1) = v(randi(length(v)));
            else % if doesn't tumble out it stays there
                V(i,1) = 0;
            end
            end
    
            T(i) = T(i)+rt+1; % resets the clock
        end
    end

    %% Calculate next positions

    Nx = Nx + V(:,1).*tstep.*cosd(A(:,2));
    Ny = Ny + V(:,1).*tstep.*sind(A(:,2));
    
    [A1,A2] = hist3([Nx,Ny],'Edges',{X_h3 Y_h3});
    
    %% Calculate next chemical
    
    g = g-A1*cons;
    
%% Extract Data
    % plot rms
    
    d = (Nx-x/2).^2 + (Ny-y/2).^2;
    rms(round(t/tstep+1),1) = sqrt(sum(d)./length(d));
    
%   refreshdata
    
    % track velocities
   
    Vrec(:,round(t/tstep+2)) = V(:,1);
   
    
    %% if you want to observe in real time
    
    % C = contourf(A2{1},A2{2},A1);
    t
    refreshdata

%     h3 = pcolor(X_h3,Y_h3,A1); %%
%     hold on  % overlays the scatter
%     h1 = scatter(Nx,Ny,25,c,'filled');
%     axis([0 x 0 y]);
%     hold off
%     h4 = pcolor(X_h3,Y_h3,g); %%
    pause(.005)
    

end
%% Analyze

    F = fit(rms,[0:tmax/tstep]','poly2');
    D(k,:) = 1./coeffvalues(F)./2;
    
    for i=1:n
        Hist_stuck(i,1) = tstep*(numel(Vrec(Vrec(i,:)==0)));
    end
    
    for t=1:tmax/tstep+2
        n_stuck(t,1) = numel(Vrec(Vrec(:,t)==0));
        n_swim(t,1) = numel(Vrec(Vrec(:,t)==60));
    end
    

        
end

%% Plots

figure()
hist(Hist_stuck,20)

trange = (0:tstep:tmax+tstep);

figure()
plot(trange,n_stuck,'r')
hold on
plot(trange,n_swim,'b')

mean(D(:,1))

%% Comments / Notes

% D can be calculated by solving this equation: <x^2> = 2*D*t
% compare this with D = vmax^2 * t / 2*(1-phi)
% 3100 vs 5400

% 5500 (stuck) vs 5400
% 
% 4930 vs 5400

% steady-state between movers and non-movers is reached in this model
