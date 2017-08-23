%% TO DO

% x. track the total number of "stuck" vs number of "swimmers"
% 2. make some stuck forever?
% x. histogram the number of times a cell is stuck
% x. incorporate "energy taxis" where cells tumble more where there is a
% high pmf (tumbling rate is higher where glucose-g is high)
% 5. incorporate a history effect on top the 
% x. calculate steadistate of stuck vs unstuck
% x. velocity decrease when in region of low food (doesn't make sense
% though, more of a gradual long-term change as the pmf is depleted)
% 6. apply reflective boundary conditions on top and bottom, apply periodic
% on sides
% 7. measure the velocity of the peak of the moving phase boundary


%% Setup
% units in micrometers
% speeds in micrometers / sec

close all
clearvars A A1 A2 c C d h1 h2 h3 Hist_stuck N0 n_stuck n_swim Nx Ny rms T V Vrec

tmax = 30; % set simulation time duration (goal is 10 mins)
tstep = 0.1; % 0.1 is standard res because the rt round goes to 0.1 resolution
iter = 1; % if repeatable

% size of area in um
x = 100;
y = 100;

D = zeros(iter,3); % creates diffusion constant measurement

for k = 1:iter
    %% housekeeping
    
    close all
    
%% Number of cells and initial condition

n = 1000; % can go up to 10e5 cells without much problem

% choose 1
% N0 = rand(n,2)*10+x/2; % center inoculation
N0 = rand(n,2)*x; % all over inoculation

% separate vectors for x and y coords
Nx = N0(:,1);
Ny = N0(:,2);

%% Velocity distribution

vmax = 20; % scales the velocity distribution
v = vmax*[1 1 1 1 1 1 1 1 1 1]; % [0 0 45 45 60 60 60 60 75 75];

for i = 1:n
    V(i,1) = v(randi(length(v)));
end

Vrec = ones(n,tmax/tstep+2);
Vrec(:,1) = V(:,1);

%% Run time distribution and rules

max_run = 3;
rt = max_run*ones(n,1);
rt0 = rt;
T = max_run*rand(n,1); % needs to be altered by local conc of glucose, glucose needs to be used up by cells
alpha = 0.6; % factor that reduces run time based on local conc

%% Tumble distribution

A = rand(n,3)*360; % the state of the tumbles
                % A(:,1) is the previous angle
                % A(:,2) is the current angle
                % A(:,3) is the previous angle it used to get stuck

athresh = 150; % angle needed to escape a dead-end

% [0 30 60 90 120 150 180];

%% Glucose field and consumption distribution

slope = 2;
cbins = 40;
g = ones(cbins,cbins);

% create a gradient of the food
for i = 1:cbins
g(i,:) = g(i,:)*2*i/cbins-1/4;
end
g(g<0) = 0;
g(g>1) = 1;
G = ones(n,1);

Cx_h4 = x/(cbins*2):x/cbins:x-x/(cbins*2);
Cy_h4 = y/(cbins*2):y/cbins:y-y/(cbins*2);

cons = 0.00005; % 0.00005

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
nbins = 40;

X_h3 = x/(nbins*2):x/nbins:x-x/(nbins*2); %A2{1}; put below
Y_h3 = y/(nbins*2):y/nbins:y-y/(nbins*2); %A2{2};

Icx = discretize(Nx,X_h3);
Icy = discretize(Ny,Y_h3);

[A1,A2] = hist3([Nx,Ny],'Edges',{X_h3 Y_h3});

% [C, h3] = contourf(A2{1},A2{2},A1);

h3 = pcolor(X_h3,Y_h3,A1);
hold on
h1 = scatter(N0(:,2),N0(:,1),25,c,'filled');
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
        
        if T(i,1)<= 0 % if the run is over
            
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
    
            T(i) = T(i)+rt(i); % resets the clock
        end
    end

%% Calculate next positions

%     V = V.*(G)+30; % V is dependent on concentration of g?
    
    Nx = Nx + V(:,1).*tstep.*cosd(A(:,2));
    Ny = Ny + V(:,1).*tstep.*sind(A(:,2));
    
    % boundary conditions
    for i=1:n
    % periodic boundary conditions on sides
    if Nx(i) < 0
        Nx(i) = x + Nx(i);
    end
    
    if Nx(i) > x
        Nx(i) = Nx(i) - x;
    end
    
    %**** UNDER CONSTRUCTION ****%
    % reflective boundary conditions on top and bottom
    if Ny(i) < 0
        Ny(i) = -Ny(i);
        if A(i,2) >= 270
        A(i,:) = A(i,:) - 2*(A(i,:)-270) - 180;
        end
        if A(i,2) < 270 && A(i,2) >= 180
        A(i,:) = A(i,:) + 2*(270-A(i,:)) - 180;
        end
    end
    
    if Ny(i) > y
        Ny(i) = Ny(i) - 2*(Ny(i) - y);
        if A(i,2) >= 90
        A(i,:) = A(i,:) - 2*(A(i,:)-90) + 180;
        end
        if A(i,2) < 90 && A(i,2) >= 0
        A(i,:) = A(i,:) + 2*(90 - A(i,:)) + 180;
        end
    end
    end
    
%% Sum number of cells per gridbox
    
    [A1,A2] = hist3([Nx,Ny],'Edges',{X_h3 Y_h3});
    
    Icx = discretize(Nx,X_h3);
    Icy = discretize(Ny,Y_h3);

    
    
%% Calculate next chemical

    g = g-A1*cons; % eating
    g(g<0) = 0;
    for i = 1:n
        if isnan(Icx(i))==0 && isnan(Icy(i))==0
            G(i,1) = g(Icx(i),Icy(i));
        end
    end
    rt = round(rt0.*(1-alpha*G),1);
    
%% Extract Data
    % plot rms
    
    d = (Nx-x/2).^2 + (Ny-y/2).^2;
    rms(round(t/tstep+1),1) = sqrt(sum(d)./length(d));
    
    % refreshdata
    
    % track velocities
   
    Vrec(:,round(t/tstep+2)) = V(:,1);
   
    
%% if you want to observe in real time
    
    % C = contourf(A2{1},A2{2},A1);
%     t
%     refreshdata

%   h3 = pcolor(X_h3,Y_h3,A1); %%
%! OFF
%     colorbar
%     hold on  % overlays the scatter % keep with below 4 on to get lines
     h1 = scatter(Ny,Nx,25,c,'filled');
     axis([0 x 0 y]);
     hold off
%     h4 = pcolor(X_h3,Y_h3,g); %%
    pause(.005)
    

end
%% Analyze

    F = fit(rms,[0:tstep:tmax]','poly2');
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
hold off

mean(D(:,1))
steady_stuck = mean(n_stuck(tmax/tstep/2:end,1))
steady_swim = n-steady_stuck

figure()
plot(rms)

figure()
pcolor(Cx_h4,Cy_h4,g)

%% Comments / Notes

% D can be calculated by solving this equation: <x^2> = 2*D*t
% compare this with D = vmax^2 * t / 2*(1-phi)
% 3100 vs 5400

% 5500 (stuck) vs 5400
% 
% 4930 vs 5400
% scaling the chance to get stuck: stuck reduces the diffusion coeff
% 0 5300 400 0
% 1 2700 235 165
% 2 1400 160 240
% 3 850 110 290
% 9 17 8 392
%
% 2^-x based decay
% exp(-x) decay

% steady-state between movers and non-movers is reached in this model
