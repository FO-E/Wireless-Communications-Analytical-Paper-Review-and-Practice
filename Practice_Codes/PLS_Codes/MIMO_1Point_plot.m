clear variables; close all; clc;

tx = 10; % Number of Transmit Antennas of MIMO
vary = 10; 
rx = linspace(1,10,vary);

%SNR_dB =linspace(1,20,vary); % SNRdB

SNR_dB = 15;
SNR = 10.^(SNR_dB/10); % SNR

N_itr = 1000;

%I = eye(N_min);


capacity_mimo = zeros(1,vary);

for i = 1 : N_itr
      
    for i = 1:vary
        
        N_min = min(i,tx); 
        I = eye(N_min);
        H = sqrt(0.5)*(randn(i,tx) + 1j*randn(i,tx));
        [~,V,~] = svd(H);
        HH = H*H';        
        capacity_mimo(i) = capacity_mimo(i) + log(real(det(I + (SNR/tx)*V(i,i)*HH )));
    
    end
      
    
end

capacity_plot = capacity_mimo/N_itr;

  
figure(1);
plot(rx,capacity_plot,'kd--');
legend('E[capacity] (N_T = 10)','location','best');
xlabel('Number of receive antennas (N_R)')
ylabel('Capacity(bps/Hz)')   
title('MIMO Capacity')
xlim([1,10]);
ylim([0,25]);
grid on;