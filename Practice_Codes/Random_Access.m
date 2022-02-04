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
%% SIMULATED RESULTS FOR PROPOSED

R = 2e3;
ep = 0.52e-6;
c = 3e8;
z = (ep*c)/2;


angles = linspace(0, 2*pi, 1000); % 1000 is the total number of points
xCenter = 0;
yCenter = 0;
x = R * cos(angles) + xCenter; 
y = R * sin(angles) + yCenter;
eNodeB = [xCenter,yCenter];
devices = 10;

prob_collision_sim = zeros(1,9);

for j = 1:9
      RA = RA_attempt(j);

      theta_1 = rand*(2*pi);
      tagged_machine_device = randsample(devices,1,true);
      r_t = sqrt(rand)*R;
      x_tagged_machine_device = tagged_machine_device + r_t.*cos(theta_1);
      y_tagged_machine_device = tagged_machine_device + r_t.*sin(theta_1);

      collision =  0;
         
         

          for n = 1:n_iteration  

                M_preambles = randperm(M);

                machines = (2*rand)/R^2;
                theta = rand(1,RA-1)*(2*pi);
                r_m = sqrt(rand(1,RA-1))*R;
                x_machine_devices = machines + r_m.*cos(theta);
                y_machine_devices = machines + r_m.*sin(theta);

                x_others = x_machine_devices;
                y_others = y_machine_devices; 

                r_o = sqrt(((x_tagged_machine_device-eNodeB(1))^2) + ((y_tagged_machine_device-eNodeB(2))^2));  %distance b/n tagged device and eNodeB

                T_o = (2*r_o)/c; % TA value for tagged device

                x_region_tagged_device = r_o * cos(angles) + xCenter; % DRAW TAGGED DEVICE CIRCLE

                y_region_tagged_device = r_o * sin(angles) + yCenter;
                

                % DRAW REGION AROUND TAGGED DEVICE

                x_region_tagged_device_1 = (r_o-eps) * cos(angles) + xCenter; 

                y_region_tagged_device_1 = (r_o-eps) * sin(angles) + yCenter;


                x_region_tagged_device_2 = (r_o+eps) * cos(angles) + xCenter; 

                y_region_tagged_device_2 = (r_o+eps) * sin(angles) + yCenter;

                tagged_device_selection = randsample(M_preambles,1,true); 
            

                     for f = 1:RA-1
                         location_comparison = 0;

                        if (sqrt((x_others(f)-eNodeB(1))^2 + (y_others(f)-eNodeB(2))^2)) >= (r_o-eps) && (sqrt((x_others(f)-eNodeB(1))^2 + (y_others(f)-eNodeB(2))^2)) <= (r_o+eps)
                           location_comparison = location_comparison + 1;
                           preamble_selection = randsample(M_preambles,1,true);
                           preamble_comparison = preamble_selection == tagged_device_selection;

                           if ((preamble_comparison == 1) && (location_comparison == 1))
                                collision = collision + 1;
%                                 break;
                          end
                        end

                     end

          end

                 prob_collision_sim(j) = collision/n_iteration;

end     

%% PLOTTING DATA POINTS

figure(1)

semilogy(k, collision_prob, 'k-', k, collision_prob_1, 'r-');
grid on; hold on; 

semilogy(k, prop_collision_prob, 'k--', k, prop_collision_prob_1, 'r--');
hold on;

semilogy(k, conv_col_prob, 'k*',k, conv_col_prob_, 'r*');
hold on;

semilogy(k, prob_collision_sim,'b*')

legend('M = 5 conv.(anal)','M = 20 conv.(anal)', 'M = 5 prop.(anal)', 'M = 20 prob.(anal)', 'M = 5 conv.(sim)', 'M = 20 conv.(sim)', 'M = 5 prop.(sim)','location', 'best')
xlabel('Number of RA attempts from machine devices') %on a single RA slot (k+1)
ylabel('Collision probability')
ylim([10^-3 10^0]);
xlim([2 10]);
