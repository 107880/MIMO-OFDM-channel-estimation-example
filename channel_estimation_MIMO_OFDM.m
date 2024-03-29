function [mse_LS, mse_LMMSE] = channel_estimation_MIMO_OFDM(Nt, Nr, H, carrier, pusch, bits, modulationScheme, snr_db)

    rep = 200; % repitition    
    subcarrier_num = 12;
    symbol_num = 14;
    pilot_num = 6;
    PRB_num = 26;
    pusch.Modulation = modulationScheme;

    % DMRS mapping indicies define
    [ind, puschInfo] = nrPUSCHIndices(carrier,pusch);
    ind = 1:subcarrier_num*PRB_num * symbol_num;
    pusch_symbols = nrPUSCH(carrier, pusch, bits);
    dmrs_tx_symbols = nrPUSCHDMRS(carrier, pusch);
    dmrs_tx_symbols_original = dmrs_tx_symbols;
    for j = 1:(Nt-1)
        dmrs_tx_symbols = [dmrs_tx_symbols; dmrs_tx_symbols_original];
    end
    
    % dmrs index processing
    dmrsInd = [];
    for sc = 1:PRB_num*subcarrier_num
        for sy = 1:symbol_num
            if(mod(sy, symbol_num) == 3)
                if( (mod(sc, subcarrier_num) ~= 13) && (mod(mod(sc, subcarrier_num), 2) == 1)  )
                    dmrsInd(end + 1) = (sy * symbol_num) + sc;
                end
            end
        end
    end
    
    
    % symbol mapping using index
    OFDM_grid = nrResourceGrid(carrier);
    OFDM_grid(ind) = pusch_symbols;
    OFDM_grid(dmrsInd) = dmrs_tx_symbols_original;
    OFDM_grid = reshape(OFDM_grid, [1, PRB_num*subcarrier_num*symbol_num]);
    OFDM_grid_flatten = repmat(OFDM_grid, [Nt, 1]);
    dmrsInd_orignal = dmrsInd;
    dmrs_tx_idx = dmrsInd_orignal;
    dmrs_rx_idx = dmrsInd_orignal;

    % dmrs index processing for tx, rx signal
    for j = 1:(Nt-1)
        dmrs_tx_idx= [dmrsInd_orignal, dmrs_tx_idx+4368];
    end
    for j = 1:(Nr-1)
        dmrs_rx_idx= [dmrsInd_orignal, dmrs_rx_idx+4368];
    end
    
    dmrs_tx_symbols = transpose(dmrs_tx_symbols);
    dmrs_tx_symbols = transpose(reshape(dmrs_tx_symbols, Nt, (pilot_num*PRB_num)));
    
    %% estimation for rayleigh channel
    % Propagation with rayleigh fading channel
    rx_signal = zeros(Nr,PRB_num*subcarrier_num*symbol_num);
    for i = 1:PRB_num*subcarrier_num*symbol_num
        X_temp = OFDM_grid_flatten(:,i);
        H_temp = squeeze(H(:,:,i));
        rx_signal(:,i) = H_temp * X_temp;
    end
    
    mse_LS = zeros(length(snr_db), 1);
    mse_LMMSE = zeros(length(snr_db), 1);
    
    for it = 1:rep
        i1 = 1;
        for snrIdx = 1:length(snr_db)
            SNR = 10^(snr_db(snrIdx) / 10);
            % no = 1.0 / SNR;
            
            noise = sqrt(1/SNR) * (randn(size(rx_signal)) + 1i*randn(size(rx_signal)));
            rx_signal_noised = (rx_signal) + noise;
            
            
            dmrs_rx_symbols = rx_signal_noised(dmrs_rx_idx);
            dmrs_rx_symbols = transpose(reshape(dmrs_rx_symbols, Nr, (pilot_num*PRB_num )));
            
    
            H_p_LS = zeros(Nr, Nt, PRB_num*pilot_num);
            H_p_LMMSE = zeros(Nr, Nt, PRB_num*pilot_num);
    
            for i=1:PRB_num*pilot_num
                X_p = dmrs_tx_symbols(i,:);
                Y_p = dmrs_rx_symbols(i,:);
                % LS estimation
                Y_inv = pinv(X_p);
                H_p_LS_temp = transpose(Y_inv * Y_p);
                H_p_LS(:,:,i) = H_p_LS_temp;
                
                % LMMSE estimation
                H_vec = reshape(H_p_LS_temp, 1, []);
                R_H = (H_vec)' * (H_vec);
                W = ( R_H / (R_H + (1/SNR) * eye(size(R_H))) );
                H_LMMSE_flatten = W * H_vec';
                H_p_LMMSE(:,:,i) = reshape(H_LMMSE_flatten, Nr, Nt);
    
            end
            % H_LS = interp1(dmrsInd, H_p_LS, 1:nfft, 'linear', 'extrap') /nrOFDMInfo(carrier).Nfft;
            % 2-dimensional interpolation
            % gridSize = [PRB_num*subcarrier_num];
            % [dmrsSubcarrier, dmrsSymbol] = ind2sub(gridSize, dmrs_rx_idx);
            % [X, Y] = meshgrid(1:gridSize(1), 1:gridSize(2));
            % H_LS = griddata(dmrsSubcarrier, dmrsSymbol, H_p_LS, X', Y', 'nearest'); 
    
        
            % calcuate MSE
            H_p = H(:,:,dmrsInd);
            diff_ls = H_p - H_p_LS;
            diff_lmmse = H_p - H_p_LMMSE;
            % diff_ls = H_rayleigh - H_LS;
            % diff_lmmse = H_rayleigh - H_LMMSE;
    
            mse_LS(i1) = mse_LS(i1) + mean(abs(diff_ls(:)) .^2);
            mse_LMMSE(i1) = mse_LMMSE(i1) + mean(abs(diff_lmmse(:)) .^2);
        
            i1 = i1+1;
        end
    end
    mse_LS = mse_LS / rep;
    mse_LMMSE = mse_LMMSE / rep;

end

