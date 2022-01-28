clear variables; close all; clc;

%% DEFINING NEEDED PARAMETERS

vary = 9;
k = linspace(2,10,vary);
M = 5;
M_1 = 20;
n_iteration = 100000;
R = 2000;
e = 0.52*10^-6;
c = 3*10^8;
eps = (e*c)/2;

%% DEFINING FUNCTIONS

fun1 = @(r) r*( 1 - ( R^2 - ((r - eps)^2) ) /(M*R^2) ).^k; 
fun2 = @(r) r*( 1 - ((r + eps)^2) /(M*R^2) ).^k; 
fun3 = @(r) r*( 1 - (4*r*eps) / (M*R^2) ).^k; 

fun_1 = @(r) r*( 1 - ( R^2 - ((r - eps)^2) ) /(M_1*R^2) ).^k; 
fun_2 = @(r) r*( 1 - ((r + eps)^2) /(M_1*R^2) ).^k; 
fun_3 = @(r) r*( 1 - (4*r*eps) / (M_1*R^2) ).^k; 

%% DATA POINTS

% collision_prob = zeros(1,vary);
% collision_prob_1 = zeros(1,vary);

%% COMPUTING PROBABILITY OF COLLISION FOR BOTH CONVENTIONAL AND PROPOSED SCHEMES

    
collision_prob = 1 - (1 - 1/M).^(k-1);
    
collision_prob_1 = 1 - (1 - 1/M_1).^(k-1);
    
prop_collision_prob = 1 - ( (2/R^2)*( integral(fun1, R-eps, R, 'ArrayValued', true) + integral(fun2, 0, eps, 'ArrayValued', true) + ....
        integral(fun3, eps, R-eps, 'ArrayValued', true) ) );
    
prop_collision_prob_1 = 1 - ( (2/R^2)*( integral(fun_1, R-eps, R, 'ArrayValued', true) + integral(fun_2, 0, eps, 'ArrayValued', true) + ....
        integral(fun_3, eps, R-eps, 'ArrayValued', true) ) );
    

%% SIMULATED RESULTS FOR CONVENTIONAL

conv_col_prob = zeros(1,9);
RA_attempt = linspace(2,10,9);

for i1=1:9
    
    RA_sample = RA_attempt(i1);
    
    preamble_matrix = randi([1, M], RA_sample, n_iteration);
    col_event = 0;
    
    for i2=1:n_iteration
        
        test = preamble_matrix(:,i2);
        
        check = sum( logical(test == test(1)) ); % Using the first column value as the tagged device
        
        if( check > 1)
            col_event = col_event + 1;
        end
              
    end
    
    conv_col_prob(i1) = col_event/n_iteration;
    
end

%% SIMULATED RESULTS FOR CONVENTIONAL

conv_col_prob_ = zeros(1,9);
RA_attempt = linspace(2,10,9);

for i1=1:9
    
    RA_sample = RA_attempt(i1);
    
    preamble_matrix = randi([1, M_1], RA_sample, n_iteration);
    col_event = 0;
    
    for i2=1:n_iteration
        
        test = preamble_matrix(:,i2);
        
        check = sum( logical(test == test(1)) ); % Using the first column value as the tagged device
        
        if( check > 1)
            col_event = col_event + 1;
        end
              
    end
    
    conv_col_prob_(i1) = col_event/n_iteration;
    
end

%% PLOTTING DATA POINTS

figure(1)

semilogy(k, collision_prob, 'k-', k, collision_prob_1, 'r-');
grid on; hold on; 

semilogy(k, prop_collision_prob, 'k--', k, prop_collision_prob_1, 'r--');
hold on;

semilogy(k, conv_col_prob, 'k*');
semilogy(k, conv_col_prob_, 'r*');

legend('M = 5 conv.(anal)','M = 20 conv.(anal)', 'M = 5 prop.(anal)', 'M = 20 prob.(anal)', 'M = 5 conv.(sim)', 'M = 20 conv.(sim)', 'location', 'best')
xlabel('Number of RA attempts from machine devices') %on a single RA slot (k+1)
ylabel('Collision probability')
ylim([10^-3 10^0]);
xlim([2 10]);
