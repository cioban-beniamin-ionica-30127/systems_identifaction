clear 
clc
load product_16.mat
plot(time,detrend(y));

m = 3;
%% Identifeicare
MSE_val = zeros(1,m);
MSE_id = zeros(1,m);
N_id = floor(4/5*length(y));
N_val = length(y) - N_id;
time_id = time(1:N_id);
y_id = detrend(y(1:N_id));
y_val = y(N_id + 1:length(y));
time_val = time(N_id + 1:length(y));

PHI = phi_finder(m,time_id);

%% Validare
PHI_val = phi_finder(m,time_val);


theta = PHI\y_id;
y_cap = PHI*theta;
y_cap_val = PHI_val*theta;
subplot(2,1,1),plot(time_id,detrend(y_id),LineWidth=2);
hold on;
subplot(2,1,1),plot(1:N_id,detrend(y_cap),'--red',LineWidth=2);

subplot(2,1,2),plot(time_val,detrend(y_val),LineWidth=2);
hold on;
subplot(2,1,2),plot(time_val,detrend(y_cap_val),'--green',LineWidth=2);

%% MSE graph

figure;
for i = 1:m

    PHI = phi_finder(i,time_id);
    PHI_val = phi_finder(i,time_val);
    theta = PHI\y_id;
    y_cap = PHI*theta;
    y_cap_val = PHI_val*theta;
    epsilon_id = y_id - y_cap;
    MSE_id(i) = 1/N_id*sum(epsilon_id.^2);
    epsilon_val = y_val - y_cap_val;
    MSE_val(i) = 1/N_val*sum(epsilon_val.^2);


end
plot(1:m,MSE_id,LineWidth=2);
hold on;
plot(1:m,MSE_val,LineWidth=2);

%% Function

function PHI = phi_finder(m,time)
    PHI = zeros(length(time),2*m+2);
    index = 1:m;
    i = 1:length(time);
    j = 3:2:2*m+1;
    PHI(i,1) = 1;
    PHI(i,2) = time;
    PHI(i,j) = cos((2*pi*index.*time)/12);
    PHI(i,j+1) = sin((2*pi*index.*time)/12);

end

