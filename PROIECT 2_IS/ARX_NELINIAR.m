clear
close all
clc
load iddata-07.mat
plot(id);

%% Generarea puterilor regresorilor (monoame de na+nb variabile)

clear powers

m = 5;
na = 1;
nb = 1;

powers = gen_powers(m,na,nb);
M = length(powers);

%% Generare matrice regresori

clear PHI y_val theta_pred
close all

PHI = get_phi(powers,id.u,id.y,na,nb);
theta_pred = PHI\id.y;
PHI_val_pred = get_phi(powers,val.u,val.y,na,nb);
%% Validare
[yhat_pred, yhat_sim] = NARX(val.u,theta_pred,na,nb,PHI_val_pred,powers);
yhat_sim = iddata(yhat_sim,val.u,val.Ts);
yhat_pred = iddata(yhat_pred,val.u,val.Ts);
figure
compare(yhat_sim,val);
figure
compare(yhat_pred,val);
mse

%% MSE
tic
clear comb
nr = 1:3;
na = 1:3;
nb = na;
m = 1:7;
N_val = length(val.y);
N_id = length(id.y);
comb = combvec(m,nr);
comb(:,[1 2]) = comb(:,[2 1]);
format longEng
mse_pred_id  = zeros(1,length(comb));
mse_sim_id = zeros(1,length(comb));
mse_pred_val = zeros(1,length(comb));
mse_sim_val = zeros(1,length(comb));
c = 1;
for i = nr
        for k = m
            format longEng
            powers = gen_powers(k,i,i);
            PHI_id_pred = get_phi(powers,id.u,id.y,i,i);
            theta_pred = PHI_id_pred\id.y;
            PHI_val_pred = get_phi(powers,val.u,val.y,i,i);
            [yhat_pred_id, yhat_sim_id] = NARX(id.u,theta_pred,i,i,PHI_id_pred,powers);
            [yhat_pred_val, yhat_sim_val] = NARX(val.u,theta_pred,i,i,PHI_val_pred,powers);
            mse_pred_id(c) = 1/N_id*sum((yhat_pred_id-id.y).^2);
            mse_sim_id(c) = 1/N_id*sum((yhat_sim_id-id.y).^2);
            mse_pred_val(c) = 1/N_val*sum((yhat_pred_val-val.y).^2);
            mse_sim_val(c) = 1/N_val*sum((yhat_sim_val-val.y).^2);
            c = c+1;
        end
end
toc

%% Best fitting and tuning

[mse_min_pred,index_pred] = min(mse_pred_val,[],'linear');
[mse_min_sim,index_sim] = min(mse_sim_val,[],'linear');
best_model_pred = comb(:,index_pred);
best_model_sim = comb(:,index_sim);

% 3D PLOT PENTRU MSE IN CAZUL PREDICTIEI
mse_pred_id_matrix = reshape(mse_pred_id,length(m),length(na));
mse_pred_val_matrix = reshape(mse_pred_val,length(m),length(na));
subplot(121), surf(na,m,mse_pred_id_matrix)
xlabel('na=nb'), ylabel('m'), zlabel('MSE')
title({'MSE','pe datele de identificare'})
subplot(122), surf(na,m,mse_pred_val_matrix)
xlabel('na=nb'), ylabel('m'), zlabel('MSE')
title({'Eroarea medie patratica','pe datele de validare'})
sgtitle('MSE pentru prediectie depinzand de na = nb si m')
%3D PLOT PENTRU MSE IN CAZUL SIMULARII
figure
mse_sim_id = rmmissing(mse_sim_id);
mse_sim_val = rmmissing(mse_sim_val);
mse_sim_id_matrix = reshape(mse_sim_id,[],length(na));
mse_sim_val_matrix = reshape(mse_sim_val,[],length(na));
subplot(121), surf(na,1:3,mse_sim_id_matrix)
xlabel('na=nb'), ylabel('m'), zlabel('MSE')
title({'MSE','pe datele de identificare'})
subplot(122), surf(na,1:3,mse_sim_val_matrix)
xlabel('na=nb'), ylabel('m'), zlabel('MSE')
title({'Eroarea medie patratica','pe datele de validare'})
sgtitle('MSE pentru simulare depinzand de na = nb si m')
%% Rezultate finale 
na_pred = best_model_pred(2); nb_pred = na_pred;

m_pred = best_model_pred(1);
na_sim = best_model_sim(2); nb_sim = na_sim;
m_sim = best_model_sim(1);

powers_pred = gen_powers(m_pred,na_pred,nb_pred);
PHI_id_pred = get_phi(powers_pred,id.u,id.y,na_pred,nb_pred);
theta_pred = PHI_id_pred\id.y;
PHI_val_pred = get_phi(powers_pred,val.u,val.y,na_pred,nb_pred);
[yhat_pred,~] = NARX(val.u,theta_pred,na_pred,nb_pred,PHI_val_pred,powers_pred);
powers_sim = gen_powers(m_sim,na_sim,nb_sim);
PHI_id_sim = get_phi(powers_sim,id.u,id.y,na_sim,nb_sim);
theta_sim = PHI_id_sim\id.y;
PHI_val_sim = get_phi(powers_sim,val.u,val.y,na_sim,nb_sim);

[~,yhat_sim] = NARX(val.u,theta_sim,na_sim,nb_sim,PHI_val_sim,powers_sim);

yhat_pred = iddata(yhat_pred,val.u,val.Ts);
compare(yhat_pred,val);
yhat_sim = iddata(yhat_sim,val.u,val.Ts);
figure
compare(yhat_sim,val);



function [PHI] = get_phi(powers,u,y,na,nb)
N = length(y);
M = length(powers);
PHI = ones(N,M);
for i = 1:N
    for j = 1 :M
        for k = 1 : na
            index = i - k ;
            if index > 0
                PHI(i,j) = PHI(i,j)*y(index)^powers(j,k);
            else
                PHI(i,j) = 0;
                break;

            end
        end

        for k = 1 : nb
            index = i - k;
            if index > 0 && PHI(i,j) ~= 0
                PHI(i,j) = PHI(i,j)*u(index)^powers(j,na+k);
            else
                PHI(i,j) = 0;
                break;
            end
        end

    end
end
end




function [powers] = gen_powers(m,na,nb)
N = na + nb;
if m > 1
    numbers = 0:m^(na+nb)-1;
    numbers = dec2base(numbers,m);
    digits_of_numbers = zeros(m^N,N);
    for i = 1 : m^N
        for M = 1:N
            digits_of_numbers(i,M) = str2double(numbers(i,M));
        end
    end
    M = 1;
    for i = 1:m^N
        if sum(digits_of_numbers(i,:)) <= m
            powers(M,:) = digits_of_numbers(i,:);
            M = M + 1;
        end
    end

    powers = [powers; m*eye(N)];
else if m == 1
        powers = eye(N);
end
end
end

function [yhat_pred, yhat_sim] = NARX(u,theta,na,nb,PHI,powers)

yhat_pred = PHI*theta;
aux = ones(length(theta),1);
[N,~] = size(PHI);
yhat_sim = zeros(N,1);
for i = 1:N
    for j = 1:length(theta)
        for k = 1 : na
            index = i - k;
            if index > 0
                aux(j) = aux(j)*yhat_sim(index)^powers(j,k);
            else
                aux(j) = 0;
                break;
            end

        end
        for k = 1 : nb
            if index > 0 && aux(j) ~= 0
                aux(j) = aux(j)*u(index)^powers(j,k+na);
            else
                aux(j) = 0;
                break;
            end
        end

    end
    yhat_sim(i) = aux'*theta;
    aux = ones(length(theta),1);
end

end

