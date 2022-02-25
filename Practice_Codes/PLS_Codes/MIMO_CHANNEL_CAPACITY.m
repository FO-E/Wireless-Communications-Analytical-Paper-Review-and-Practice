clear variables; close all; clc;

rx = 4; % Number of Receive Antennas of MIMO 
tx = 20; % Number of Transmit Antennas of MIMO
tx_new = 10;

N_min = min(rx,tx); 
N_min_new = min(rx,tx_new);

vary = 20; 
SNR_dB =linspace(1,20,vary); % SNRdB
SNR = 10.^(SNR_dB/10); % SNR
N_itr = 1000;

I = eye(N_min);
I_new = eye(N_min_new);

capacity_mimo = zeros(1,vary);
capacity_mimo_new = zeros(1,vary);

for i = 1 : N_itr
    
    H = sqrt(0.5)*(randn(rx,tx) + 1j*randn(rx,tx));
    H_new = sqrt(0.5)*(randn(rx,tx_new) + 1j*randn(rx,tx_new));
    
    [~,V,~] = svd(H);
    [~,V_new,~] = svd(H_new);
    
    HH = H*H';
    H_newH_new = H_new*H_new';
  
    for i = 1:vary
                
        capacity_mimo(i) = capacity_mimo(i) + log(real(det(I + (SNR(i)/tx)*V(2,2)*HH )));
        capacity_mimo_new(i) = capacity_mimo_new(i) + log(real(det(I_new + (SNR(i)/tx_new)*V_new(3,3)*H_newH_new )));
    
    end
      
    
end

capacity_plot = capacity_mimo/N_itr;
capacity_plot_new = capacity_mimo_new/N_itr;
  
figure(1);
plot(SNR_dB,capacity_plot,'kx--', SNR_dB,capacity_plot_new,'r*--');
legend('E[capacity] (N_T = 20)','E[capacity] (N_T = 10)','location','best');
xlabel('SNR (dB)')
ylabel('Capacity(bps/Hz)')   
title('MIMO Capacity')
xlim([5,20]);
ylim([0,30]);
grid on;