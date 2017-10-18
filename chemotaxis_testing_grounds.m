%% Single Cell Chemotaxis Testing Grounds

tic

% no consumption
% http://matlabtricks.com/post-44/generate-random-numbers-with-a-given-distribution

%% MAIN CONTROLS

record_vid = 0; % Set to 1 if you want to record a video
t_framegrab = 5;

rt_graph_setting = 'Agents';
cells_show = 1; % Set to 1 if you want to show individual cells
chemo1 = 0; % Energy taxis !!!
chemo2 = 1; % Classic chemotaxis !!!

labels = {'density_3d' 'n' 'x' 'y' 'tmax' 'tstep' 'vmax' 'pstuck' 'athresh' 'max_run' 'min_run' 'alpha' 'cons' 'slope' 'avg_peak_vel (um/s)' 'avg_peak_vel (mm/hr)' 'R2' 'avg_peak_size' 'winsize'};

folder = 'chemotaxis_testing_grounds/';
name = 'test';
mkdir(strcat(folder,name));

%% Control Graph Output

iter = 1; % length(stuck_mat); % repeatable

color_conc = 0; % 1 to turn on to bin conc density and graph it
max_A1 = 12; % max cell number in grid box

% 'RMS' if you want to see the root mean square displacement

rotational_diffusion = 1;

% Set to 1 if you want rotational diffusion
% Set to 0 if you want to save computational time

%% Setup
% units in micrometers
% speeds in micrometers / sec

tmax = 500; % set simulation time duration (goal is 10 mins) %$ %100
tstep = 0.05; % 0.1 is standard res because the rt round goes to 0.1 resolution %$

% size of area in um
x = 400; % 400 %$
y = 4000; % 4000 %$ %800

if strcmp(rt_graph_setting,'RMS')
D = zeros(iter,3); % creates diffusion constant measurement
end

%% Housekeeping
    
    close all
    clearvars A A1 A1_t A2 c C d h1 h2 h3 Hist_stuck N0 n_stuck n_swim Nx Ny rms peaks T V Vrec
    
%% Number of Cells and Initial Condition

density_3d = 0.001; % set the 3d density of cells (assuming slice size of 5 um) % 0.0025 or 4 times that

n = density_3d*x*y*5 % can go up to 10e5 cells without much problem % 10000 % 100

% choose 1 layout of cells

% 1 % center inoculation
% N0 = rand(n,2)*10+x/2; 

% 2 % evenly distributed all over space, merely append *x
N0 = rand(n,2);
N0(:,1) = N0(:,1)*x;
N0(:,2) = N0(:,2)*y;

% separate vectors for x and y coords
Nx = N0(:,1);
Ny = N0(:,2);

%% Velocity Distribution

vmax = 20; % scales the velocity distribution % 20 %$
pstuck = 0 ; % 0.8 %$
v = vmax*[1 1 1 1 1 1 1 1 1 1 1]; % [0 0 45 45 60 60 60 60 75 75];
istuck = round(length(v)*pstuck); % appends velocity distribution with stuck
v(1,end+1:end+istuck) = zeros(1,istuck);

V = zeros(n,2);
for i = 1:n
    V(i,1) = v(randi(length(v))); % column 1 can be zero
end
V(:,2) = V(:,1); % column 2 saves the original velocity

%% Run Time Distribution and Rules

% see chemo_response.m
% https://www.desmos.com/calculator/tvb2scjlbw

ave_run = 0.6; %$
rt = round(exprnd(ave_run,n,1),1); % next run time, changes based on g
T = round(exprnd(ave_run,n,1),1); % initial tumble count-down vector
convolution = zeros(n,1); % initial convolution state !!!

%% Tumble Angle Distribution

A = rand(n,3)*360;
if rotational_diffusion
rot_dif = zeros(n,1);
D_rot = 1; % Hz %$
end

%% Glucose Field and Consumption Distribution

slope = 150; %100 %$ % y/25 y/5 %30
cbins.x = x; %20 %$ %x/20
cbins.y = y; %400 %$ %y/20
g = ones(cbins.y,cbins.x);

% create a gradient of the food across y
for i = 1:cbins.y
g(i,:) = g(i,:)*slope*i/cbins.y; % -50
end
g(g<0) = 0;
g(g>slope*y/2) = slope*y/2;

if chemo1 %!!!
G = ones(n,1)*slope*y/2;
end

if chemo2
cwin = ceil(10*ave_run); %!!!
G = zeros(n,cwin/tstep); %!!!
end

Cx_h4 = x/(cbins.x*2):x/cbins.x:x-x/(cbins.x*2);
Cy_h4 = y/(cbins.y*2):y/cbins.y:y-y/(cbins.y*2);

%% SIMULATION

%% Initialize Real-time Graphs

% figure() % initialize rms evolution
if strcmp(rt_graph_setting,'RMS')
rms = zeros(tmax/tstep+1,1);
h2 = plot(0:tstep:tmax,rms);
axis([0,tmax,0,1000])
h2.YDataSource = 'rms';
end

c = [rand(n,1) rand(n,1) rand(n,1)]; % random colors for each agent
nbins.x = cbins.x;
nbins.y = cbins.y;

X_h3 = 0:x/nbins.x:x; %A2{1}; put below
Y_h3 = 0:y/nbins.y:y; %A2{2};

Icx = discretize(Ny,X_h3); % switched because matrix format
Icy = discretize(Nx,Y_h3);

%% Tumbling Angle Distribution

%!!!
theta = -170:10:180;
%data = ?
%p_theta = interp1(length(data),data,theta,'spline')
pdf_theta = 1/2*(1+cosd(theta)); % psuedomonas: 1/16*(3+15*cosd(theta)*cosd(theta));
pdf_theta = pdf_theta / sum(pdf_theta);

cdf_theta = cumsum(pdf_theta);
[cdf_theta, mask] = unique(cdf_theta);
theta = theta(mask);
rng_tum = interp1(cdf_theta,theta,rand(1,25));

for t=0:tstep:tmax
    
%% Construct State Vectors (stuck, moving, tumbling)
 
    T= T-tstep;

    %idx = T<=0;
    
    for i=1:n            
        if T(i,1)<= 0 % if the run is over            
            rng_tum = interp1(cdf_theta,theta,rand(1)); % input the experimental angle distribution
            
            A(i,2) = wrapTo360(A(i,2) + rng_tum); % conditional change in direction
            
            T(i) = T(i)+rt(i); % resets the clock
        end
        
        if isnan(T(i))
            error('T has received an NaN')
        end
    end

%% Calculate Next Positions

    % rotational diffusion
if rotational_diffusion
    rot_dif = normrnd(0,1,n,1)*sqrt(2*D_rot*tstep);
end

    Nx = Nx + V(:,1).*tstep.*cosd(A(:,2)+rot_dif);
    Ny = Ny + V(:,1).*tstep.*sind(A(:,2)+rot_dif);
    
    % boundary conditions
    for i=1:n
    % periodic boundary conditions on sides
    if Nx(i) < 0
        Nx(i) = x + Nx(i);
    end
    
    if Nx(i) > x
        Nx(i) = Nx(i) - x;
    end
    
    
    % reflective boundary conditions on top and bottom
    if Ny(i) < 0
        Ny(i) = -Ny(i);
        if A(i,2) >= 270 && A(i,2) < 360
        A(i,:) = wrapTo360(A(i,:) - 2*(A(i,:)-270) - 180);
        end
        if A(i,2) <= 270 && A(i,2) > 180
        A(i,:) = wrapTo360(A(i,:) + 2*(270-A(i,:)) - 180);
        end
    end
    
    if Ny(i) > y
        Ny(i) = Ny(i) - 2*(Ny(i) - y);
        if A(i,2) >= 90 && A(i,2) < 180
        A(i,:) = wrapTo360(A(i,:) - 2*(A(i,:)-90) + 180);
        end
        if A(i,2) <= 90 && A(i,2) > 0
        A(i,:) = wrapTo360(A(i,:) + 2*(90-A(i,:)) + 180);
        end
    end
    end
    
%% Calculate Next Chemical
    
if chemo2
    G = circshift(G,[0 -1]); %!!!
end

Icx = discretize(Ny,Y_h3,'IncludedEdge','right'); % switched because matrix format
Icy = discretize(Nx,X_h3,'IncludedEdge','right');

    for i = 1:n % in case any cells fall out of the sim they have a G to use
        if isnan(Icx(i))==0 && isnan(Icy(i))==0
            G(i,end) = g(Icx(i),Icy(i));
        end
    end
    
%% Chemotaxis

    min_run = 0.2; %$
    max_run = 5;

if chemo2
    convolution = zeros(n,1); % reset convolution state !!!
    for t1=tstep:tstep:cwin %!!!
        convolution = convolution + tstep*chemo_response(cwin-t1,ave_run,0.018).*G(:,round(t1/tstep));
    end
    ave_run_chem = ave_run./(1-convolution); %!!! 
    ave_run_chem(ave_run_chem>max_run) = max_run;
    ave_run_chem(ave_run_chem<min_run) = min_run;
    rt = exprnd(ave_run_chem); %!!!
    rt(rt>max_run) = max_run;
    rt(rt<min_run) = min_run;
end
    
if chemo1
    rt0 = round(exprnd(ave_run,n,1),1);
    rt = round(rt0.*(1-alpha*G),1); % tumbling rate response to chemical
    rt(rt>max_run) = max_run;
    rt(rt<min_run) = min_run;
end
    
    
%% Extract Data
    % plot rms in real time
if strcmp(rt_graph_setting,'RMS')
    d = (Nx-x/2).^2 + (Ny-y/2).^2;
    rms(round(t/tstep+1),1) = sqrt(sum(d)./length(d));
    refreshdata
end   
    
%% Update Real-Time Graphs for the Cells in Real Time

if mod(t,t_framegrab)==0
if strcmp(rt_graph_setting,'Agents')
    
    if color_conc
    h4 = pcolor(X_h3,Y_h3,g); %% turn off to not update concentration    
    colormap bone
    colorbar
    caxis([0 1])
    end
    
     if cells_show
     %! hold on  % overlays the scatter % keep with below 4 on to get lines
     h1 = scatter(Nx,Ny,25,c,'filled');
     %! hold off
     end

     title(strcat(num2str(t),'sec'))
     axis([0 x 0 y]);
     pbaspect([1 y/x 1]);
     pause(0.001)
end
end
    
if record_vid
    frame(round(t/t_framegrab+1)) = getframe(gcf);
end

end

toc