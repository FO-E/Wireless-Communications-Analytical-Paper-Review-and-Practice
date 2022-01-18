clear variables; close all; clc;

%% DEFINING NEEDED PARAMETERS

vary = 9;
k = linspace(2,10,vary);
M = 5;
M_1 = 20;
N_Itr = 100;
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

collision_prob = zeros(1,vary);
collision_prob_1 = zeros(1,vary);

prop_collision_prob = zeros(1,vary);
prop_collision_prob_1 = zeros(1,vary);

%% COMPUTING PROBABILITY OF COLLISION FOR BOTH CONVENTIONAL AND PROPOSED SCHEMES

for i = 1:1:N_Itr
    
    prob_col = 1 - (1 - 1/M).^(k-1);
    
    prob_col_1 = 1 - (1 - 1/M_1).^(k-1);
    
    prop_prob_col = 1 - ( (2/R^2)*( integral(fun1, R-eps, R, 'ArrayValued', true) + integral(fun2, 0, eps, 'ArrayValued', true) + ....
        integral(fun3, eps, R-eps, 'ArrayValued', true) ) );
    
    prop_prob_col_1 = 1 - ( (2/R^2)*( integral(fun_1, R-eps, R, 'ArrayValued', true) + integral(fun_2, 0, eps, 'ArrayValued', true) + ....
        integral(fun_3, eps, R-eps, 'ArrayValued', true) ) );
    
    collision_prob = collision_prob + prob_col;
    
    collision_prob_1 = collision_prob_1 + prob_col_1;
    
    prop_collision_prob = prop_collision_prob + prop_prob_col; 
    
    prop_collision_prob_1 = prop_collision_prob_1 + prop_prob_col_1;
    
end

collision_prob = collision_prob / N_Itr;
collision_prob_1 = collision_prob_1 / N_Itr;
prop_collision_prob = prop_collision_prob / N_Itr;
prop_collision_prob_1 = prop_collision_prob_1 / N_Itr;

%% PLOTTING DATA POINTS

figure(1)

semilogy(k, collision_prob, 'k-', k, collision_prob_1, 'r-');
% legend('M = 5 conv.(anal)','M = 20 conv.(anal)', 'location' ,'best')
grid on; hold on; 
semilogy(k, prop_collision_prob, 'k--', k, prop_collision_prob_1, 'r--');
legend('M = 5 conv.(anal)','M = 20 conv.(anal)', 'M = 5 prop.(anal)', 'M = 20 prob.(anal)', 'location', 'best')

xlabel('Number of RA attempts from machine devices') %on a single RA slot (k+1)
ylabel('Collision probability')
ylim([10^-3 10^0]);
xlim([2 10]);


