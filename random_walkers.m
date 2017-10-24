%% IMPORTANT INFO TO READ BEFORE YOU RUN THIS CODE

% Create a folder in MATLAB called active_walkers_output and that will store the data
% that you tell the program to use and the results (figures, csv) it
% outputs. The tune-able parameters are marked with %$, so to find them
% simply hit ctrl-F and search %$.

% Everything is in units of microns and seconds.

%% MAIN CONTROLS

record_vid = 0; % Set to 1 if you want to record a video
t_framegrab = 10; % nth number of sim frame to record

rt_graph_setting = 'MF'; %#

% Set to 'MF' if you want to see the mean field
% 'Agents' if you want to see the agents

cells_show = 0; % Set to 1 if you want to show individual cells

labels = {'density_3d' 'n' 'x' 'y' 'tmax' 'tstep' 'vmax' 'pstuck' 'athresh' 'max_run' 'min_run' 'alpha' 'cons' 'slope' 'avg_peak_vel (um/s)' 'avg_peak_vel (mm/hr)' 'R2' 'avg_peak_size' 'winsize'};

folder = 'active_walkers_output/';
name = 'test_1';
mkdir(strcat(folder,name));

file = strcat(folder,name,'/',rt_graph_setting,'_',name,'.csv'); % creates writeable CSV
fid = fopen(file, 'w') ;
 fprintf(fid, '%s,', labels{1,1:end-1}) ;
 fprintf(fid, '%s\n', labels{1,end}) ;
 fclose(fid) ;

%% TO DO / DONE LIST

% x. track the total number of "stuck" vs number of "swimmers"
% x. make some stuck forever?
%       as a natural consequence yes
%       incorporated nonetheless
% x. histogram the number of times a cell is stuck
% x. incorporate "energy taxis" where cells tumble more where there is a
% high pmf (tumbling rate is higher where glucose-g is high)
% 2. incorporate a long-term history effect on energy taxis
% x. chemotactic response function that alters tumbling frequency
%   have finite history of chemical encountered by each cell
% for chemotaxis cut off the 
% x. calculate steadistate of stuck vs unstuck
% x. velocity decrease when in region of low food (doesn't make sense
% though, more of a gradual long-term change as the pmf is depleted)
% x. apply reflective boundary conditions on top and bottom, apply periodic
% on sides
% x. measure the velocity of the peak of the moving phase boundary
% (superimpose the cross-sections)
% x. apply the longer boundary case (4000 um in length by 200 um)
% x. how much does being stuck influence the velocity of this phase
% boundary?
%       not much
% 5. measure stuck vs swimming time and space dependence in and near band
% given the constraints, how long do they get stuck for?
% x. implement easy settings for switching off graphs or components
% x. get a csv file as an output
% x. organize photos and csv's into a folder
% x. match the velocity distributions to data (append 0 velocities at end)
% 7. create an input matrix that determines initial conditions of sim
% 8. make it so that bacteria that travel faster have a higher chance of
% getting stuck
% 9. power law on tumbling frequency
% x. control density instead of cell number
% 10. incorporate autochemotaxis and see what it does
% x. diffusion on chemicals
% 11. save last band configuration (Nx Ny convolution rt)
% 12. bin dependence?
% 13. calculate CMC, value of ave_run_chem in space and time
%       https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3164075/

%% Real World 

% majority of cells are not motile 5-25%
% fit data to pdf
% power law or gamma distribution tumbling rate
% poisson run velocity
% cells maintain run velocity across runs
% random tumble angle
% chemotactic response function
% glucose chemotaxis
% consumption and diffusion of chemicals
% higher run speed gives more directional persistence = https://www.nature.com/articles/ncomms9776#supplementary-information


%% Control Graph Output

iter = 1; % length(stuck_mat); % repeatable

color_cell = 1; % 1 to turn on to bin cell density and graph it
color_conc = 0; % 1 to turn on to bin conc density and graph it
max_A1 = 12; % max cell number in grid box

% 'RMS' if you want to see the root mean square displacement

end_graph_stuck = 0;

% Set to 1 if you want to see the histogram of how long cells are stuck
% Set to 0 if you want to save memory

output_steady = 1;

% Set to 1 if you want to know how many cells are swimming or not
% Set to 0 to hide output

rotational_diffusion = 1;

% Set to 1 if you want rotational diffusion
% Set to 0 if you want to save computational time

velocity_dep = 0;

% Set to 1 if you want velocity to depend on environment
% Set to 0 if you want predetermined velocity

batch_mode = 0;

% Set to 1 if you want to iterate a bunch of parameters in a matrix
% Set to 0 to manually do things one at a time

% Choose tumbling probability modulation parameter
chemo1 = 0; % Energy taxis
chemo2 = 1; % Classic chemotaxis

%% Setup
% units in micrometers
% speeds in micrometers / sec

tmax = 260; % set simulation time duration (goal is 10 mins) %$ %100
tstep = 0.1; % 0.1 is standard res because the rt round goes to 0.1 resolution %$

% size of area in um
x = 800; % 400 %$
y = 4000; % 4000 %$

if strcmp(rt_graph_setting,'RMS')
D = zeros(iter,3); % creates diffusion constant measurement
end

for k = 1:iter
%% Housekeeping
    
    close all
    clearvars A A1 A1_t A2 c C d h1 h2 h3 Hist_stuck N0 n_stuck n_swim Nx Ny rms peaks T V Vrec
    
%% Number of Cells and Initial Condition

density_3d = 0.0001; % set the 3d density of cells (assuming slice size of 5 um) % 0.0025 or 4 times that

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

peaks = zeros(tmax/tstep+1,2); % recording matrix with peak value and loci

%% Velocity Distribution

pdf_v = [];
v_ave = 20;
v = v_ave/10:v_ave/10:v_ave*2;
pdf_v = normpdf(v,v_ave,v_ave/3);
% vmax = 20; % scales the velocity distribution % 20 %$
immobile = 0; % 0.1 0.8 %$  % percentage of cells that are immobile at beginning
%v = vmax*[1 1 1 1 1 1 1 1 1 1 1]; % [0 0 45 45 60 60 60 60 75 75];

if immobile > 0
pdf_v = [immobile/(1-immobile) pdf_v]; % appends velocity distribution with stuck
v = [0 v];
end

pdf_v = pdf_v / sum(pdf_v);
cdf_v = cumsum(pdf_v);
[cdf_v, mask] = unique(cdf_v);
v = v(mask);

V = zeros(n,2);
V(:,1) = interp1(cdf_v,v,rand(n,1),'nearest'); % column 1 can be zero
V(isnan(V)) = 0;
V(:,2) = V(:,1); % column 2 saves the original velocity

% probability to get stuck
stuck = 0; % see data %0.2
if stuck > 0
Vpick = [zeros(n,10*stuck) V(:,2).*ones(n,10*(1-stuck))];
end

if end_graph_stuck
Vrec = ones(n,tmax/tstep+2);
Vrec(:,1) = V(:,1);
end

%% Run Time Distribution and Rules

% see chemo_response.m
% https://www.desmos.com/calculator/tvb2scjlbw

if chemo1
ave_run = 0.8; %$
alpha = 0.6; % the larger this number, the larger the factor that reduces run time based on local conc % 0.6 %$
end

if chemo2
ave_run = 0.6;
convolution = zeros(n,1); % initial convolution state for chemotactic response function
end

rt = round(exprnd(ave_run,n,1),1); % next run time, changes based on g
% rt0 = rt;
T = round(exprnd(ave_run,n,1),1); % initial tumble count-down vector

%% Tumble Angle Distribution

athresh = 150; % angle needed to escape a dead-end %$

A = rand(n,3)*360;
                % A(:,1) is the previous absolute angle at previous time
                % A(:,2) is the current absolute angle
                % A(:,3) is the previous angle it used to get stuck

if rotational_diffusion
rot_dif = zeros(n,1);
D_rot = 1; % Hz %$
end

theta = -170:10:180;
%data = ?
%p_theta = interp1(length(data),data,theta,'spline')
pdf_theta = 1/2*(1+cosd(theta)); % psuedomonas: 1/16*(3+15*cosd(theta)*cosd(theta));
pdf_theta = pdf_theta / sum(pdf_theta);

cdf_theta = cumsum(pdf_theta);
[cdf_theta, mask] = unique(cdf_theta);
theta = theta(mask);
theta = [-180 theta];
cdf_theta = [0 cdf_theta];
%rng_tum = interp1(cdf_theta,theta,rand(1,25));

%% Glucose Field and Consumption Distribution

Dg = 100; % Diffusion %200

binsize = 20; %2
slope = 0.05; %100 %$ % y/25 y/5 %30 %0.1 %0.02
cbins.x = x/binsize; %20 %$ %x/20
cbins.y = y/binsize; %400 %$ %y/20
g = ones(cbins.y,cbins.x);
g_max = slope*y;

% create a gradient of the food across y %!
for i = 1:cbins.y
g(i,:) = g(i,:)*slope*i*y/cbins.y - g_max/4; % -50 %!!!
end
g(g<0) = 0;
g(g>g_max/4) = g_max/4;

if chemo1 %!!!
G = ones(n,1)*slope*y/2;
end

if chemo2
cwin = ceil(10*ave_run); %!!!
G = zeros(n,cwin/tstep); %!!!
end

Cx_h4 = 0:x/cbins.x:x;
Cy_h4 = 0:y/cbins.y:y;

cons = 0; % 0.00005 %$ %0.005 %0.002

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

X_h3 = 0:x/nbins.x:x-x/nbins.x; %A2{1}; put below
Y_h3 = 0:y/nbins.y:y-y/nbins.y; %A2{2};

    Icx = discretize(Ny,[Y_h3 y]); % switched because matrix format
    Icy = discretize(Nx,[X_h3 x]);

[A1,A2] = hist3([Nx,Ny],'Edges',{X_h3 Y_h3});

% SCRAPS [C, h3] = contourf(A2{1},A2{2},A1);

if strcmp(rt_graph_setting,'Agents')

figure() % scatter and color plot

if cells_show == 0
h3 = pcolor(X_h3,Y_h3,A1'); % A1 is transposed to match x y coords
else
h4 = pcolor(X_h3,Y_h3,g); %% turn off to not update concentration
hold on % superimpose agents on top of field
h1 = scatter(N0(:,1),N0(:,2),25,c,'filled');
axis([0 x 0 y]);
hold off
end

h1.XDataSource = 'Nx';
h1.YDataSource = 'Ny';

end

if strcmp(rt_graph_setting,'MF') 
figure() % mean field output initalize
h5 = plot(mean(g));
hold on
h6 = plot(mean(A1));
%title(strcat(num2str(t),'sec'))
xlabel('Distance (um)')
ylabel('Cell Density (orange) & Food Density (blue)')
hold off
end

for t=0:tstep:tmax
    
%% Construct State Vectors (stuck, moving, tumbling)
 
        % countdown
        % velocity
        % angle
    T= T-tstep;
    
    
    
    for i=1:n
            
        % A(:,1) = A(:,2);
        %%%%%% Get indices should be faster
        
        if T(i,1)<= 0 % if the run is over
            
            rng_tum = interp1(cdf_theta,theta,rand(1),'nearest'); % input the experimental angle distribution    
            
            A(i,1) = A(i,2); % sets the history
            
            if V(i,1) ~= 0
                A(i,3) = A(i,2); % if cell is in the stuck state, it will maintain its A(i,3) value
            end
            
            A(i,2) = wrapTo360(A(i,2) + rng_tum); % conditional change in direction
            
            % agar tunnel condition
            if V(i,1) ~= 0
            V(i,1) = Vpick(i,randi(length(Vpick(1,:)))); % pick random velocity in distribution

            else % must have a certain angle to escape
            if abs(A(i,2)-A(i,3)) >= athresh && abs(A(i,2)-A(i,3)) <= athresh+(180-athresh)*2 && V(i,2) ~= 0
                V(i,1) = V(i,2);
            end
            end
    
            T(i) = T(i)+rt(i); % resets the clock
        
        % Errors
            
        if isnan(T(i)) %!!!
            error('T has received an NaN')
        end
        
        if isnan(A(i,2)) %!!!
            error('A has received an NaN')
        end
        
        end
    end

%% Calculate Next Positions

if velocity_dep
     V = V.*(G)+30; % V is dependent on concentration of g? incorporate time delay
end

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
    
%% Sum up the Number of Cells per Gridbox
    
    [A1,A2] = hist3([Nx,Ny],'Edges',{X_h3 Y_h3});
    
    Icx = discretize(Ny,[Y_h3 y]); % switched because matrix format
    Icy = discretize(Nx,[X_h3 x]);
    
%% Calculate Next Chemical

    g = g-A1'*cons; % eating
    g(g<0) = 0;
    
   %!!!
    for i = 1:n % in case any cells fall out of the sim they have a G to use
        if isnan(Icx(i))==0 && isnan(Icy(i))==0
            G(i,1) = g(Icx(i),Icy(i));
        end
    end
    
    % Diffusion
    
    % Shifted stuff for x
    
    sgr = circshift(g, [0 1]);
    sgl = circshift(g, [0 -1]);
    sgd = circshift(g, [1 0]);
    sgu = circshift(g, [-1 0]);
    
    % Laplacian term for x
    
    dx = 20;
    
    laplgx = (sgr+sgl-2*g)/(dx)/(dx);     
    laplgy = (sgd+sgu-2*g)/(dx)/(dx);
    laplg = laplgx + laplgy;
    
    g = g + Dg*laplg*tstep;
    g(end,:) = g(end-1,:);
    g(1,:) = g(2,:);
    
%% Chemotaxis
    
    min_run = 0; %$ %0.2
    max_run = 5; %$
    
if chemo1
    rt0 = round(exprnd(ave_run,n,1),1);
    rt = round(rt0.*(1-alpha*G),1); % tumbling rate response to chemical
    rt(rt>max_run) = max_run;
    rt(rt<min_run) = min_run;
end

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
    
%% Extract Data
    % plot rms in real time
    if strcmp(rt_graph_setting,'RMS')
    d = (Nx-x/2).^2 + (Ny-y/2).^2;
    rms(round(t/tstep+1),1) = sqrt(sum(d)./length(d));
    refreshdata
    end

    % track velocities to get distribution of stuck vs non-stuck
    if end_graph_stuck
    Vrec(:,round(t/tstep+2)) = V(:,1);
    end
   
    
%% Update Real-Time Graphs for the Cells in Real Time

A1_t = mean(A1)';
winsize = 10; % prepare to filter the noise with a moving window average %$
b = (1/winsize)*ones(1,winsize);
a = 1;
A1_filt = filter(b,a,A1_t);

[peak, peak_loci] = max(A1_filt); % extracts the peak
peaks(round(t/tstep+1),1) = peak_loci*y/cbins.y;
peaks(round(t/tstep+1),2) = peak;

% mean field output to 1D cross-section
if mod(t,t_framegrab)==0
if strcmp(rt_graph_setting,'MF')
h5 = plot(Y_h3,mean(g(:,1:end-1)'));
title(strcat(num2str(t),'sec'))
ylim([0 n/cbins.x/cbins.y*2])
xlim([0 y])
xlabel('Distance (um)')
ylabel('Cell Density (orange) & Food Density (blue)')

hold on
h6 = plot(Y_h3,A1_filt); % A1_filt if you want to see the smoother version
h7 = scatter(peak_loci*y/cbins.y, peak,50,'r','filled');
hold off

pause(0.001)
end
end

if mod(t,t_framegrab)==0
if strcmp(rt_graph_setting,'Agents')
    if color_cell
    h3 = pcolor(X_h3,Y_h3,A1'); %% turn off to not update the cell density color
    colormap bone
    colorbar
    caxis([0 max_A1])
    end
    if color_conc
    h4 = pcolor(X_h3,Y_h3,g); %% turn off to not update concentration    
    colormap bone
    colorbar
    caxis([0 1])
    end
  
     if cells_show
     hold on  % overlays the scatter % keep with below 4 on to get lines
     h1 = scatter(Nx,Ny,25,c,'filled');
     hold off
     end

     title(strcat(num2str(t),'sec'))
     axis([0 x 0 y]);
     pbaspect([1 y/x 1]);
     ylabel('Microns')
     xlabel('Microns')
     pause(0.001)
end
end
    
if record_vid
    if mod(t,t_framegrab)==0
    frame(round(t/t_framegrab+1)) = getframe(gcf);
    end
end

end
%% Analyze Data

if strcmp(rt_graph_setting,'RMS')
    % Diffusion coefficient fit
F = fit(rms,[0:tstep:tmax]','poly2');
D(k,:) = 1./coeffvalues(F)./2;
end

% Find the portion that is stuck verses swimming
if end_graph_stuck
    for i=1:n
        Hist_stuck(i,1) = tstep*(numel(Vrec(Vrec(i,:)==0)));
    end
    
    for t=1:tmax/tstep+2
        n_stuck(t,1) = numel(Vrec(Vrec(:,t)==0));
        n_swim(t,1) = numel(Vrec(Vrec(:,t)==20));
    end
end    


%% Plots

% how long on average a cell is stuck for?
if end_graph_stuck
figure()
hist(Hist_stuck,20)
trange = (0:tstep:tmax+tstep);

% is a steady state reached, and how much?
figure()
plot(trange,n_stuck,'r')
hold on
plot(trange,n_swim,'b')
hold off

% tabulates portion of stuck cells at steady state
steady_stuck = mean(n_stuck(tmax/tstep/2:end,1));
steady_swim = n-steady_stuck;
stuck_frac = steady_stuck / n;
end

mint = 200;

if strcmp(rt_graph_setting,'RMS')
mean(D(:,1))
figure()
plot(rms)
end

m = polyfit(mint:tmax/tstep, peaks(mint:tmax/tstep), 1); % find average velocity
R = corrcoef(mint:tmax/tstep, peaks(mint:tmax/tstep));
R2 = R.^2;
R2 = R2(2,1); % R-squared value

if color_conc
figure()
pcolor(Cx_h4,Cy_h4,g)
title('Concentration at Final Time')
colorbar
pbaspect([1 y/x 1]);
end

avg_peak_vel = m(1,1);
avg_peak_size = mean(peaks(mint:end,2)); % average height

figure()
plot(peaks(:,1))
title('Peak Motion over Time')
saveas(gcf,strcat(folder,name,'/',rt_graph_setting,'_',name,'_',num2str(avg_peak_vel),'_peak.jpg'))


%% Write CSV Report File

%! UNDER CONSTRUCTION

labels = {'density_3d' 'n' 'x' 'y' 'tmax' 'tstep' 'v_ave' 'pstuck' 'athresh' 'max_run' 'min_run' 'alpha' 'cons' 'slope' 'avg_peak_vel (um/s)' 'avg_peak_vel (mm/hr)' 'R2' 'avg_peak_size' 'winsize'};
values = [density_3d n x y tmax tstep v_ave pstuck athresh max_run min_run alpha cons slope avg_peak_vel avg_peak_vel*3.6 R2 avg_peak_size winsize];
dlmwrite(file, values, '-append') ; % write params

if end_graph_stuck
stuck_frac
end

if strcmp(rt_graph_setting,'MF')
avg_peak_vel
avg_peak_size
end
        
end

if record_vid
    vidfile = strcat(folder,'/',name,'/video_',name);
    codec = 'MPEG-4';
    video = VideoWriter(vidfile,codec);
    open(video)
    writeVideo(video,frame)
    close(video)
end

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

% a steady-state between movers and non-movers is reached in this model
% 0.5 to 1.0 um/s band speed is the goal
% 
% cells occupy about 0.004 of space in microscope
% 0.0015 cells per um^2
%

%% Beta 0.1 Results

% changing alpha does not change band speed appreciably, but raising alpha
% raises avg peak height slightly
% 
% decreasing overall frequency of tumbles resulted in higher peaks, but did
% not affect band speed as much
%
% slope of gradient seems to be the greatest determining factor of speed of
% peak by e^-slope
%
% Cells adjust tumbling rate depending on level of food in environment

%% Beta 0.2 Results

% Changes: chemotaxis response function
% Cells adjust tumbling rate depending on gradient of food in environment

%% Beta 0.3 Results

% 