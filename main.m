% Title: MIMO-OFDM channel estimation simulator
% date: 2024.03.29
% name: 107880
% description: conventional channel estimator in MIMO-OFDM

clear all
clc

carrierFrequency = 4e9;
rmsDelaySpread = rand * (300e-9 - 10e-9) + 10e-9;
maxDopplerShift = rand * (500 - 0) + 0; 
snrDb = rand * (32 + 4) - 4;
sirDb = rand * (36 - 0) + 0;
subcarrierSpacing = 15; % [GHz]
ofdmSymbolDuration = 71e-6;
ttiLength = 1;
codeRate = 658/1024;
nfft = 128;
path_num = 16;

Nt = 4;
Nr = 8;
subcarrier_num = 12;
symbol_num = 14;
pilot_num = 6;
PRB_num = 26;
modulationScheme = '16QAM';
symbol_bit = 4;



%% TX 
% generate signal and modulation
bits = randi([0 1], PRB_num * subcarrier_num * symbol_num * symbol_bit, 1);

bits_QPSK = randi([0 1], PRB_num * subcarrier_num * symbol_num * 2, 1);
bits_16QAM = randi([0 1], PRB_num * subcarrier_num * symbol_num * 4, 1);
bits_64QAM = randi([0 1], PRB_num * subcarrier_num * symbol_num * 6, 1);


% OFDM parameter: carrier config
carrier = nrCarrierConfig;
carrier.SubcarrierSpacing = subcarrierSpacing;
carrier.CyclicPrefix = 'normal';
carrier.NSizeGrid = PRB_num;
OFDMinfo = nrOFDMInfo(carrier);
% carrier.SymbolsPerSlot = 14; % read only


% PUSCH setting
p = 1; % # of layers
pusch = nrPUSCHConfig('NumLayers',p);
pusch.PRBSet = 0:(PRB_num-1);



%% channel modeling
% rayleigh fading channel
H_rayleigh = (randn(Nr, Nt, PRB_num*subcarrier_num*symbol_num) + 1i*randn(Nr, Nt, PRB_num*subcarrier_num*symbol_num)) / sqrt(2);
H_rayleigh_1x1 = (randn(2, 2, PRB_num*subcarrier_num*symbol_num) + 1i*randn(2, 2, PRB_num*subcarrier_num*symbol_num)) / sqrt(2);
H_rayleigh_4x4 = (randn(4, 4, PRB_num*subcarrier_num*symbol_num) + 1i*randn(4, 4, PRB_num*subcarrier_num*symbol_num)) / sqrt(2);
H_rayleigh_16x16 = (randn(16, 16, PRB_num*subcarrier_num*symbol_num) + 1i*randn(16, 16, PRB_num*subcarrier_num*symbol_num)) / sqrt(2);


% rician fading channel(LOS)
K = 3; % LOS / NLOS
LOS = sqrt(K/(K+1));
NLOS = (randn(Nr, Nt, PRB_num*subcarrier_num*symbol_num) + 1i*(randn(Nr, Nt, PRB_num*subcarrier_num*symbol_num))) / sqrt(2*(K+1));
H_rician = NLOS + LOS;


% exponential correlation
rho = 0.5;
Rt = rho.^(abs((1:Nt) - (1:Nt).'));
H_exp = zeros(Nr, Nt, PRB_num*subcarrier_num*symbol_num);
for k = 1:PRB_num*subcarrier_num*symbol_num;
    H_w = (randn(Nr, Nt) + 1i*randn(Nr, Nt)) / sqrt(2);
    H_exp(:,:,k) = H_w * sqrtm(Rt);
end



%% channel estimation
snr_db = -10:5:20;
[mse_LS_rayleigh, mse_LMMSE_rayleigh] = channel_estimation_MIMO_OFDM(Nt, Nr, H_rayleigh, carrier, pusch, bits, modulationScheme, snr_db);
[mse_LS_rician, mse_LMMSE_rician] = channel_estimation_MIMO_OFDM(Nt, Nr, H_rician, carrier, pusch, bits, modulationScheme, snr_db);
[mse_LS_exp, mse_LMMSE_exp] = channel_estimation_MIMO_OFDM(Nt, Nr, H_exp, carrier, pusch, bits, modulationScheme, snr_db);

[mse_LS_QPSK, mse_LMMSE_QPSK] = channel_estimation_MIMO_OFDM(Nt, Nr, H_rayleigh, carrier, pusch, bits_QPSK, 'QPSK', snr_db);
[mse_LS_16QAM, mse_LMMSE_16QAM] = channel_estimation_MIMO_OFDM(Nt, Nr, H_rayleigh, carrier, pusch, bits_16QAM, '16QAM', snr_db);
[mse_LS_64QAM, mse_LMMSE_64QAM] = channel_estimation_MIMO_OFDM(Nt, Nr, H_rayleigh, carrier, pusch, bits_64QAM, '64QAM', snr_db);

[mse_LS_2, mse_LMMSE_2] = channel_estimation_MIMO_OFDM(2, 2, H_rayleigh_1x1, carrier, pusch, bits, modulationScheme, snr_db);
[mse_LS_4, mse_LMMSE_4] = channel_estimation_MIMO_OFDM(4, 4, H_rayleigh_4x4, carrier, pusch, bits, modulationScheme, snr_db);
[mse_LS_16, mse_LMMSE_16] = channel_estimation_MIMO_OFDM(16, 16, H_rayleigh_16x16, carrier, pusch, bits, modulationScheme, snr_db);


    
%% plot result: comparision of Rayleigh, Rician, Exponential correlation channel

max_list = [mse_LMMSE_exp, mse_LS_exp, mse_LMMSE_rician, mse_LS_rician, mse_LMMSE_rayleigh, mse_LS_rayleigh];
max_value = max(max(max_list));

semilogy(snr_db, mse_LS_rayleigh/(max_value), '-or', 'LineWidth',2);
hold on;
semilogy(snr_db, mse_LMMSE_rayleigh/(max_value), '-.or', 'LineWidth',2);
semilogy(snr_db, mse_LS_rician/(max_value), '-vb', 'LineWidth',2);
semilogy(snr_db, mse_LMMSE_rician/(max_value), '-.vb', 'LineWidth',2);
semilogy(snr_db, mse_LS_exp/(max_value), '-diamondk', 'LineWidth',2);
semilogy(snr_db, mse_LMMSE_exp/(max_value), '-.diamondk', 'LineWidth',2);

grid on;
title(['Nr = ' num2str(Nr), ' Nt = ' num2str(Nt), ': MIMO-OFDM channel estimation']);
xlabel('SNR[dB]'); 
ylabel('Normalized Channel MSE');
legend('Rayleigh-LS', 'Rayleigh-LMMSE', 'Rician-LS', 'Rician-LMMSE', 'ExpCorr-LS', 'ExpCorr-LMMSE');


%% plot result: comparision of modulation: QPSK, 16QAM, 64QAM

max_list = [mse_LMMSE_QPSK, mse_LS_QPSK, mse_LMMSE_16QAM, mse_LS_16QAM, mse_LMMSE_64QAM, mse_LS_64QAM];
max_value = max(max(max_list));

semilogy(snr_db, mse_LS_QPSK/(max_value), '-or', 'LineWidth',2);
hold on;
semilogy(snr_db, mse_LMMSE_QPSK/(max_value), '-.or', 'LineWidth',2);
semilogy(snr_db, mse_LS_16QAM/(max_value), '-vb', 'LineWidth',2);
semilogy(snr_db, mse_LMMSE_16QAM/(max_value), '-.vb', 'LineWidth',2);
semilogy(snr_db, mse_LS_64QAM/(max_value), '-diamondk', 'LineWidth',2);
semilogy(snr_db, mse_LMMSE_64QAM/(max_value), '-.diamondk', 'LineWidth',2);

grid on;
title(['Nr = ' num2str(Nr), ' Nt = ' num2str(Nt), ': MIMO-OFDM channel estimation']);
xlabel('SNR[dB]'); 
ylabel('Normalized Channel MSE');
legend('QPSK-LS', 'QPSK-LMMSE', '16QAM-LS', '16QAM-LMMSE', '64QAM-LS', '64QAM-LMMSE');


%% plot result: comparision of # of antenna(Nr, Nt)

max_list = [mse_LMMSE_2, mse_LS_2, mse_LMMSE_4, mse_LS_4, mse_LMMSE_16, mse_LS_16];
max_value = max(max(max_list));

semilogy(snr_db, mse_LS_2/(max_value), '-or', 'LineWidth',2);
hold on;
semilogy(snr_db, mse_LMMSE_2/(max_value), '-.or', 'LineWidth',2);
semilogy(snr_db, mse_LS_4/(max_value), '-vb', 'LineWidth',2);
semilogy(snr_db, mse_LMMSE_4/(max_value), '-.vb', 'LineWidth',2);
semilogy(snr_db, mse_LS_16/(max_value), '-diamondk', 'LineWidth',2);
semilogy(snr_db, mse_LMMSE_16/(max_value), '-.diamondk', 'LineWidth',2);

grid on;
title(['MIMO-OFDM channel estimation: ', modulationScheme]);
xlabel('SNR[dB]'); 
ylabel('Normalized Channel MSE');
legend('2x2-LS', '2x2-LMMSE', '4x4-LS', '4x4-LMMSE', '16x16-LS', '16x16-LMMSE');




