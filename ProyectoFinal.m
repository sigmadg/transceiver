 fileTx = 'upiita.jpg';
 txGain = -20;

[rxData, txData, bandera] = imagen(fileTx, txGain)

    % Realizar el cálculo de la tasa de error de bits (BER)
    bitErrorRate = comm.ErrorRate;
    err = bitErrorRate(double(rxData(:)), txData(1:length(reshape(rxData,[],1))));
    fprintf('  \n Tasa de error de bit (BER):\n');
    fprintf('  Tasa de error de bit (BER) = %0.5f.\n',err(1));
    fprintf('  Número de errores en los bits = %d.\n', err(2));
    fprintf('  Número de bits transmitidos = %d.\n\n',length(txData));
    
    if bandera == 2
        if err(2) >= 1
            fprintf(' El programa tiene muchos errores \n\n');
            [rxData, txData, bandera] = imagen_1(fileTx, txGain)
            
            % Realizar el cálculo de la tasa de error de bits (BER)
            bitErrorRate = comm.ErrorRate;
            err = bitErrorRate(double(rxData(:)), txData(1:length(reshape(rxData,[],1))));
            fprintf(' \n Tasa de error de bit (BER):\n');
            fprintf('  Tasa de error de bit (BER) = %0.5f.\n',err(1));
            fprintf('  Número de errores en los bits = %d.\n', err(2));
            fprintf('  Número de bits transmitidos = %d.\n\n',length(txData));
            
             if err(2) >= 1
                fprintf(' El programa tiene muchos errores \n\n');
                [rxData, txData, bandera] = ProyectoFinal(fileTx, txGain)
             else
                fprintf(' Ok \n\n');
             end
            
        else
            fprintf(' Ok \n\n');
        end
    else
        if err(2) >= 1
            fprintf(' El programa tiene muchos errores \n\n');
            [rxData, txData, bandera] = ProyectoFinal(fileTx, txGain)
        else
            fprintf(' Ok \n\n');
        end
   end