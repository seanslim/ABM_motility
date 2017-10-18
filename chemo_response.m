% t is integrated time, ave_run is run length in absence of
% chemoattractants
% https://www.desmos.com/calculator/tvb2scjlbw
function R = chemo_response(t,ave_run,a)
t0 = 1.25; %1.2
y0 = 0.023824; %0.39112
beta = 9; %9.4
R = a/ave_run.*(1-beta*((t-t0)/ave_run + (t-t0).^2/(2*ave_run.^2))).*exp((-t+t0)/ave_run) - y0;
end

%% Test for zero integration
% test = chemo_response(0:0.1:20,0.6,1)
% sum(test)
% plot(test)
