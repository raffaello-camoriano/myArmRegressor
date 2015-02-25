clear all;
close all;

%% Load ground truth
datasetFileName = 'part1-left.csv';
dataset = importdata(datasetFileName, ',', 3);

skip = 0;
f_gnd = dataset.data(skip+1:end,13:15);
t_gnd = dataset.data(skip+1:end,16:18);

n =  size(f_gnd,1);

%% Load predictions
predFileName = 'results.csv';
pred = load(predFileName);

f_est = pred(:,1:3);
t_est = pred(:,4:6);

f_cad = pred(:,7:9);
t_cad = pred(:,10:12);


%% Compute averaged RMSEs for cad and estimated parametric predictions at each step

% Estimated parameters

avg_mse_f_est = zeros(n,1);
avg_mse_t_est = zeros(n,1);
avg_mse_f_cad = zeros(n,1);
avg_mse_t_cad = zeros(n,1);

for i = 1:n
    % RIVEDERE!!!   
    if i > 1
        avg_mse_f_est(i) = (i-1)/i * avg_mse_f_est(i-1) + mean((f_est(i,:) - f_gnd(i,:)).^2 / i , 2);
        avg_mse_t_est(i) = (i-1)/i * avg_mse_t_est(i-1) + mean((t_est(i,:) - t_gnd(i,:)).^2 / i , 2);
        avg_mse_f_cad(i) = (i-1)/i * avg_mse_f_cad(i-1) + mean((f_cad(i,:) - f_gnd(i,:)).^2 / i , 2);
        avg_mse_t_cad(i) = (i-1)/i * avg_mse_t_cad(i-1) + mean((t_cad(i,:) - t_gnd(i,:)).^2 / i , 2);
    else
        avg_mse_f_est(i) = mean((f_est(i,:) - f_gnd(i,:)).^2 , 2);
        avg_mse_t_est(i) = mean((t_est(i,:) - t_gnd(i,:)).^2 , 2);        
        avg_mse_f_cad(i) = mean((f_cad(i,:) - f_gnd(i,:)).^2 , 2);
        avg_mse_t_cad(i) = mean((t_cad(i,:) - t_gnd(i,:)).^2 , 2);        
    end    
end

avg_rmse_f_est = sqrt(avg_mse_f_est);
avg_rmse_t_est = sqrt(avg_mse_t_est);
avg_rmse_f_cad = sqrt(avg_mse_f_cad);
avg_rmse_t_cad = sqrt(avg_mse_t_cad);

%% Plot

figure
title('Output 1')
hold on
plot (f_gnd(:,1))
plot (f_est(:,1))
plot (f_cad(:,1))
legend('Ground truth','Estimated pars.', 'CAD pars.')
hold off

figure
title('Output 2')
hold on
plot (f_gnd(:,2))
plot (f_est(:,2))
plot (f_cad(:,2))
legend('Ground truth','Estimated pars.', 'CAD pars.')
hold off

figure
title('Output 3')
hold on
plot (f_gnd(:,3))
plot (f_est(:,3))
plot (f_cad(:,3))
legend('Ground truth','Estimated pars.', 'CAD pars.')
hold off

figure
title('Output 4')
hold on
plot (t_gnd(:,1))
plot (t_est(:,1))
plot (t_cad(:,1))
legend('Ground truth','Estimated pars.', 'CAD pars.')
hold off

figure
title('Output 5')
hold on
plot (t_gnd(:,2))
plot (t_est(:,2))
plot (t_cad(:,2))
legend('Ground truth','Estimated pars.', 'CAD pars.')
hold off

figure
title('Output 6')
hold on
plot (t_gnd(:,3))
plot (t_est(:,3))
plot (t_cad(:,3))
legend('Ground truth','Estimated pars.', 'CAD pars.')
hold off

figure
hold on
plot (avg_rmse_f_est)
plot (avg_rmse_f_cad)
legend('avg rmse f est','avg rmse f cad')
hold off

figure
hold on
plot (avg_rmse_t_est)
plot (avg_rmse_t_cad)
legend('avg rmse t est','avg rmse t cad')
hold off

%savemultfigs
