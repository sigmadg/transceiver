function [rxData, txData, bandera] = imagen_1(fileTx, txGain)


if isempty(ver('wlan')) % Compruebe la instalación de WLAN Toolbox
    error('Instale WLAN Toolbox para ejecutar este ejemplo.');
elseif ~license('test', 'WLAN_System_Toolbox') % Verifique que haya una licencia válida
    error('Se requiere una licencia válida para WLAN Toolbox para ejecutar este ejemplo.');
end

% Configura todos los parametros y figuras que se mostrarán a lo largo del ejemplo.

% Controlador de configuración para el trazado de la imagen
if ~exist('imFig', 'var') || ~ishandle(imFig)
    imFig = figure;
    imFig.NumberTitle = 'off';
    imFig.Name = 'Graficar imagen';
    imFig.Visible = 'off';
else
    clf(imFig); % Limpiar figura
    imFig.Visible = 'off';
end

% fileTx = 'upiita.jpg';
% txGain = -20;

fData = imread(fileTx); % Leer datos de imagen del archivo

% Grafica la imagen transmitida
figure(imFig);
imFig.Visible = 'on';
subplot(211);
    imshow(fData);
    title('Imagen Transmitida');
subplot(212);
    title('La imagen recibida aparecerá aquí...');
    set(gca,'Visible','off');
    set(findall(gca, 'type', 'text'), 'visible', 'on');
pause(1); % Pausa para trazar la imagen transmitida


% Configurar la ventana del espectro
spectrumScope = dsp.SpectrumAnalyzer('SpectrumType', 'Power density', 'SpectralAverages', 10, 'YLimits', [-130 -50], 'Title', 'Espectro de señal WLAN de banda base recibida', 'YLabel', 'Densidad espectral de potencia','Position', [69 376 800 450]);

% Configurar el visor de diagramas de constelaciones para símbolos WLAN ecualizados
constellation = comm.ConstellationDiagram('Title', 'Símbolos de WLAN ecualizados', 'ShowReferenceConstellation', false, 'Position', [878 376 460 460]);

% El objeto del sistema <matlab:plutoradiodoc('commsdrtxpluto') se utiliza con PlutoSDR para transmitir datos de banda base al hardware SDR.

% Inicializar dispositivo SDR
deviceNameSDR = 'Pluto'; % Establecer dispositivo SDR
radio = sdrdev(deviceNameSDR); % Crea el objeto del dispositivo SDR
bandera = 1;
% txGain = -20;

% Preparar el archivo de imagen

% Ingresa un archivo de imagen y lo convierte en flujo binario
fileTx = 'upiita.jpg';             % Nombre del archivo de imagen
fData = imread(fileTx);            % Leer datos de imagen del archivo
scale = 0.2;                       % Factor de escala de imagen
origSize = size(fData);            % Tamaño de la imagen de entrada original
scaledSize = max(floor(scale.*origSize(1:2)),1); % Calcular el nuevo tamaño de la imagen
heightIx = min(round(((1:scaledSize(1))-0.5)./scale+0.5),origSize(1));
widthIx = min(round(((1:scaledSize(2))-0.5)./scale+0.5),origSize(2));
fData = fData(heightIx,widthIx,:); % Cambiar el tamaño de la imagen
imsize = size(fData);              % Almacenar un nuevo tamaño de imagen
binData = dec2bin(fData(:),8);     % Convierta a 8 bits unsigned binario
txImage = reshape((binData-'0').',1,[]).'; % Crear flujo binario

% Cuando el archivo de imagen se recibe y decodifica correctamente, muestra la imagen.

% Grafica la imagen transmitida
figure(imFig);
imFig.Visible = 'on';
subplot(211);
    imshow(fData);
    title('Imagen Transmitida');
subplot(212);
    title('La imagen recibida aparecerá aquí...');
    set(gca,'Visible','off');
    set(findall(gca, 'type', 'text'), 'visible', 'on');
pause(1); % Pausa para trazar la imagen transmitida

% El archivo de imagen | txImage |.
% El flujo de datos binarios se divide en unidades de transmisión de tamaño | msduLength |.
% Se agregan un encabezado MAC y un campo CRC a cada unidad de transmisión para constituir una MPDU.
% La longitud del MPDU no debe exceder los 4095 bytes.

msduLength = 4048; % Longitud de MSDU en bytes
msduBits = msduLength*8;
numMSDUs = ceil(length(txImage)/msduBits);
padZeros = msduBits-mod(length(txImage),msduBits);
txData = [txImage; zeros(padZeros,1)];

% El FCS se calcula utilizando el estándar polinomio generador de grado 32 como se define en la sección 8.2.4.8 del estándar 802.11n HT

generatorPolynomial = [32 26 23 22 16 12 11 10 8 7 5 4 2 1 0];
fcsGenerator = comm.CRCGenerator(generatorPolynomial);
fcsGenerator.InitialConditions = 1;
fcsGenerator.DirectMethod = true;
fcsGenerator.FinalXOR = 1;

% Divida el flujo de datos de entrada en fragmentos

numFragment = 0;
bitsPerOctet = 8;
lengthMACheader = 256;                           % Longitud del encabezado MPDU en bits
lengthFCS = 32;                                  % Longitud de FCS en bits
lengthMPDU = lengthMACheader+msduBits+lengthFCS; % Longitud de MPDU en bits
data = zeros(lengthMPDU*numMSDUs,1);

for ind=0:numMSDUs-1
   % Extrae bits para cada MPDU
   frameBody = txData(ind*msduBits+1:msduBits*(ind+1),:);
   % Generar bits de encabezado MPDU
   mpduHeader = helperNonHTMACHeader(mod(numFragment, 16),mod(ind,4096));
   % Cree MPDU con encabezado, cuerpo y FCS
   psdu = fcsGenerator([mpduHeader;frameBody]);
   % Concatenar PSDU para generar formas de onda
   data(lengthMPDU*ind+1:lengthMPDU*(ind+1)) = psdu;
end

%  Generar señal WLAN de banda base IEEE 802.11a

% Las propiedades del objeto contiene la configuración.
% En este ejemplo, un objeto es configurado para un ancho de banda de 20 MHz, 1 antena de transmisión y velocidad 64QAM 2/3 (MCS 6).

nonHTcfg = wlanNonHTConfig;                    % Crear configuración de paquetes
nonHTcfg.MCS = 6;                              % Modulación: 64QAM Tasa: 2/3
nonHTcfg.NumTransmitAntennas = 1;              % Número de antena transmisora
chanBW = nonHTcfg.ChannelBandwidth;
nonHTcfg.PSDULength = lengthMPDU/bitsPerOctet; % Establecer la longitud de PSDU

% El sdrTransmitter utiliza el comando | transmitRepeat | para transmitir la forma de onda WLAN de banda base en un bucle desde la memoria DDR en PlutoSDR.
% La señal de RF transmitida se sobremuestrea y se transmite a 30 MHz.
% La señal% 802.11a se transmite en el canal 5, que corresponde a un centro de frecuencia de 2.432 GHz según se define en la sección 17.4.6.3 de [1].

sdrTransmitter = sdrtx(deviceNameSDR); % Propiedades del transmisor
sdrTransmitter.RadioID = 'usb:0';

% Vuelva a muestrear la forma de onda de transmisión a 30 MHz

fs = wlanSampleRate(nonHTcfg); % Transmitir frecuencia de muestreo en MHz
osf = 1.5;                     % Factor de sobremuestreo

sdrTransmitter.BasebandSampleRate = fs*osf;
sdrTransmitter.CenterFrequency = 1e9;  % Canal 5
sdrTransmitter.ShowAdvancedProperties = true;
sdrTransmitter.Gain = txGain;

% Inicializa el codificador con un número entero aleatorio para cada paquete
scramblerInitialization = randi([1 127],numMSDUs,1);

% Genere paquetes NonHT de banda base separados por tiempo de inactividad
txWaveform = wlanWaveformGenerator(data,nonHTcfg, 'NumPackets',numMSDUs,'IdleTime',20e-6, 'ScramblerInitialization',scramblerInitialization);

% Remuestrear la forma de onda de transmisión
txWaveform  = resample(txWaveform,fs*osf,fs);

fprintf('\nGenerando forma de onda de transmisión WLAN:\n')

% Escale la señal normalizada para evitar la saturación de las etapas de RF
powerScaleFactor = 0.8;
txWaveform = txWaveform.*(1/max(abs(txWaveform))*powerScaleFactor);

% Transmitir forma de onda de RF
sdrTransmitter.transmitRepeat(txWaveform);

% La función | transmitRepeat | transfiere los paquetes WLAN de banda base con de tiempo de inactividad al PlutoSDR y almacena las muestras de señal en el hardware de memoria.

% Configuración del receptor

% La frecuencia de muestreo del receptor es de 30 MHz, que es 1,5 veces la frecuencia de muestreo de banda base de 20 MHz.
% El objeto del sistema <matlab: plutoradiodoc ('commsdrrxpluto') SDR Receiver> es utilizado con PlutoSDR para recibir datos de banda base desde el hardware SDR.

sdrReceiver = sdrrx(deviceNameSDR);
sdrReceiver.RadioID = 'usb:0';
sdrReceiver.BasebandSampleRate = sdrTransmitter.BasebandSampleRate;
sdrReceiver.CenterFrequency = sdrTransmitter.CenterFrequency;
sdrReceiver.GainSource = 'Manual';
sdrReceiver.Gain = 10;
sdrReceiver.OutputDataType = 'double';

% Configuración para recibir muestras equivalentes al doble de la longitud de la señal transmitida, esto es para asegurar que las PSDU se reciban en orden.
% En la recepción, se eliminan los fragmentos MAC duplicados.

samplesPerFrame = length(txWaveform);
sdrReceiver.SamplesPerFrame = samplesPerFrame*2;
spectrumScope.SampleRate = sdrReceiver.BasebandSampleRate;

% Obtenga los índices de campo requeridos dentro de una PSDU

indLSTF = wlanFieldIndices(nonHTcfg,'L-STF');
indLLTF = wlanFieldIndices(nonHTcfg,'L-LTF');
indLSIG = wlanFieldIndices(nonHTcfg,'L-SIG');
Ns = indLSIG(2)-indLSIG(1)+1; % Número de muestras en un símbolo OFDM

% Capturar paquetes de recepción

% La forma de onda transmitida se captura utilizando PlutoSDR.
fprintf('\nIniciar una nueva captura de RF.\n')
burstCaptures = sdrReceiver();

% Procesamiento del receptor

% Se usa un bucle while para capturar y decodificar paquetes.
% Cuando el paquete WLAN se ha decodificado correctamente, muestra el número de la secuencia detectada en la ventana de comandos para cada paquete recibido.

% Muestra la densidad espectral de potencia de la forma de onda recibida
spectrumScope(burstCaptures);

% Reducir la muestra de la señal recibida
rxWaveform = resample(burstCaptures,fs,fs*osf);
rxWaveformLen = size(rxWaveform,1);
searchOffset = 0; % Desplazamiento desde el inicio de la forma de onda en muestras

% La longitud mínima del paquete es de 10 símbolos OFDM
lstfLen = double(indLSTF(2)); % Número de muestras en L-STF
minPktLen = lstfLen*5;
pktInd = 1;
sr = wlanSampleRate(nonHTcfg); % Tasa de muestreo
fineTimingOffset = [];
packetSeq = [];
displayFlag = 0; % Bandera para mostrar la información decodificada

% Genera FCS para MPDU
fcsDetector = comm.CRCDetector(generatorPolynomial);
fcsDetector.InitialConditions = 1;
fcsDetector.DirectMethod = true;
fcsDetector.FinalXOR = 1;

% Realizar el cálculo de EVM
evmCalculator = comm.EVM('AveragingDimensions',[1 2 3]);
evmCalculator.MaximumEVMOutputPort = true;

% Procesamiento del receptor

while (searchOffset + minPktLen) <= rxWaveformLen
    % Detección de paquetes
    pktOffset = wlanPacketDetect(rxWaveform, chanBW, searchOffset, 0.8);

    % Ajustar la compensación de paquetes
    pktOffset = searchOffset+pktOffset;
    if isempty(pktOffset) || (pktOffset+double(indLSIG(2))>rxWaveformLen)
        if pktInd==1
            disp('** No se detectó ningún paquete **');
        end
        break;
    end

    % Extraiga campos que no sean HT y realice una corrección de desplazamiento de frecuencia aproximada para permitir una sincronización de símbolos confiable.
    nonHT = rxWaveform(pktOffset+(indLSTF(1):indLSIG(2)),:);
    coarseFreqOffset = wlanCoarseCFOEstimate(nonHT,chanBW);
    nonHT = helperFrequencyOffset(nonHT,fs,-coarseFreqOffset);

    % Sincronización de tiempo de símbolo
    fineTimingOffset = wlanSymbolTimingEstimate(nonHT,chanBW);

    % Ajustar la compensación de paquetes
    pktOffset = pktOffset+fineTimingOffset;

    % Extrae el campo de preámbulo NonHT después de la sincronización y realiza corrección de frecuencia
    if (pktOffset<0) || ((pktOffset+minPktLen)>rxWaveformLen)
        searchOffset = pktOffset+1.5*lstfLen;
        continue;
    end
    fprintf('\nPaquete-%d detectado en el índice %d\n',pktInd,pktOffset+1);

    % Extrai los primeros 7 símbolos OFDM para la detección y decodificación L-SIG
    nonHT = rxWaveform(pktOffset+(1:7*Ns),:);
    nonHT = helperFrequencyOffset(nonHT,fs,-coarseFreqOffset);

    % Realiza la corrección de frecuencia en la sincronización
    lltf = nonHT(indLLTF(1):indLLTF(2),:);           % Extraer L-LTF
    fineFreqOffset = wlanFineCFOEstimate(lltf,chanBW);
    nonHT = helperFrequencyOffset(nonHT,fs,-fineFreqOffset);
    cfoCorrection = coarseFreqOffset+fineFreqOffset; % Total CFO

    % Estimación de canal usando L-LTF
    lltf = nonHT(indLLTF(1):indLLTF(2),:);
    demodLLTF = wlanLLTFDemodulate(lltf,chanBW);
    chanEstLLTF = wlanLLTFChannelEstimate(demodLLTF,chanBW);

    % Estimación de ruido
    noiseVarNonHT = helperNoiseEstimate(demodLLTF);

    % Detección de formato de paquete usando los 3 símbolos OFDM inmediatamente después de L-LTF
    format = wlanFormatDetect(nonHT(indLLTF(2)+(1:3*Ns),:), chanEstLLTF,noiseVarNonHT,chanBW);
    disp(['  ' format ' formato detectado']);
    if ~strcmp(format,'Non-HT')
        fprintf('  Se ha detectado un formato diferente a Non-HT\n');
        searchOffset = pktOffset+1.5*lstfLen;
        continue;
    end

    % Recuperar bits de campo L-SIG
    [recLSIGBits,failCheck] = wlanLSIGRecover(nonHT(indLSIG(1):indLSIG(2),:), chanEstLLTF,noiseVarNonHT,chanBW);

    if failCheck
        fprintf('  Fallo de verificación L-SIG \n');
        searchOffset = pktOffset+1.5*lstfLen;
        continue;
    else
        fprintf('  Aprobacion de verificación L-SIG \n');
    end

    % Recuperar parámetros de paquetes basados en L-SIG decodificado
    [lsigMCS,lsigLen,rxSamples] = helperInterpretLSIG(recLSIGBits,sr);

    if (rxSamples+pktOffset)>length(rxWaveform)
        disp('** No hay suficientes muestras para decodificar el paquete **');
        break;
    end

    % Aplicar la corrección de CFO a todo el paquete
    rxWaveform(pktOffset+(1:rxSamples),:) = helperFrequencyOffset(rxWaveform(pktOffset+(1:rxSamples),:),fs,-cfoCorrection);

    % Crear un objeto de configuración de recepción no HT
    rxNonHTcfg = wlanNonHTConfig;
    rxNonHTcfg.MCS = lsigMCS;
    rxNonHTcfg.PSDULength = lsigLen;

    % Obtiene los índices de campo de datos dentro de la PPDU
    indNonHTData = wlanFieldIndices(rxNonHTcfg,'NonHT-Data');

    % Recupere bits de PSDU utilizando parámetros de paquetes transmitidos y estimaciones de canal de L-LTF
    [rxPSDU,eqSym] = wlanNonHTDataRecover(rxWaveform(pktOffset+(indNonHTData(1):indNonHTData(2)),:), chanEstLLTF,noiseVarNonHT,rxNonHTcfg);

    constellation(reshape(eqSym,[],1)); %Grafica la constelación
    pause(0);                           %Vuelve a graficar la constelación
    release(constellation);             %Liberar la trama de la constelación

    refSym = helperReferenceSymbols(eqSym,rxNonHTcfg);
    [evm.RMS,evm.Peak] = evmCalculator(refSym,eqSym);

    % Eliminar FCS del encabezado MAC y el cuerpo de la trama
    [rxBit{pktInd},crcCheck] = fcsDetector(double(rxPSDU)); %#ok<SAGROW>

    if ~crcCheck
         disp(' Verificación de MAC CRC');
    else
         disp('  Error de verificación MAC CRC');
    end

    % Recibe fragmentos de MAC, y recupera la información de secuenciación de los fragmentos de MAC
    [mac,packetSeq(pktInd)] = helperNonHTMACHeaderDecode(rxBit{pktInd}); %#ok<SAGROW>

    % Muestra la información decodificada
    if displayFlag
        fprintf('  CFO estimado: %5.1f Hz\n\n',cfoCorrection); %#ok<UNRCH>

        disp('  Contenido de L-SIG decodificado: ');
        fprintf('                            MCS: %d\n',lsigMCS);
        fprintf('                         Longitud: %d\n',lsigLen);
        fprintf('    Número de muestras en el paquete: %d\n\n',rxSamples);

        fprintf('  EVM:\n');
        fprintf('    EVM pico: %0.3f%%  EVM RMS: %0.3f%%\n\n', ...
        evm.Peak,evm.RMS);

        fprintf('  Contenido del campo de la secuencia MAC decodificada:\n');
        fprintf('    Secuencia de números:%d\n',packetSeq(pktInd));
    end

    % Actualizar índice de búsqueda
    searchOffset = pktOffset+double(indNonHTData(2));

    pktInd = pktInd+1;
    % Finaliza el procesamiento cuando se detecta un paquete duplicado.
    % Los datos recuperados incluyen bits de tramas duplicadas
    if length(unique(packetSeq))<length(packetSeq)
        break
    end
end

% Libera el estado del objeto sdrTransmitter y sdrReceiver
release(sdrTransmitter);
release(sdrReceiver);

% *Reconstruct Image*

% La imagen se reconstruye a partir de las tramas MAC recibidas.
if ~(isempty(fineTimingOffset)||isempty(pktOffset))&& (numMSDUs==(numel(packetSeq)-1))
    % Elimina el encabezado MAC y duplica el fragmento MAC capturado
    rxBitMatrix = cell2mat(rxBit);
    rxData = rxBitMatrix(lengthMACheader+1:end,1:numel(packetSeq)-1);

    startSeq = find(packetSeq==0);
    rxData = circshift(rxData,[0 -(startSeq(1)-1)]); % Solicitar fragmentos de MAC

    % Realizar el cálculo de la tasa de error de bits (BER)
    bitErrorRate = comm.ErrorRate;
    err = bitErrorRate(double(rxData(:)), txData(1:length(reshape(rxData,[],1))));
    fprintf('  \n Tasa de error de bit (BER):\n');
    fprintf('          Tasa de error de bit (BER) = %0.5f.\n',err(1));
    fprintf('          Número de errores en los bits = %d.\n', err(2));
    fprintf('  Número de bits transmitidos = %d.\n\n',length(txData));

    % Recreate image from received data
    fprintf('\nConstruyendo una imagen a partir de los datos recibidos.\n');

    str = reshape(sprintf('%d',rxData(1:length(txImage))),8,[]).';
    decdata = uint8(bin2dec(str));

    receivedImage = reshape(decdata,imsize);
    % Grafica imagen recibida
    if exist('imFig', 'var') && ishandle(imFig) %Verifica si la figura de Tx está abierta
        figure(imFig); subplot(212);
    else
        figure; subplot(212);
    end
    imshow(receivedImage);
    title(sprintf('Imagen Recibida'));

    if err(2) >= 10
      fprintf(' El programa tiene muchos errores \n\n');

    else
      fprintf(' Ok \n\n');
    end
end

displayEndOfDemoMessage(mfilename)

end
