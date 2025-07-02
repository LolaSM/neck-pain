%% 
% *ACCESO A LAS CARPETAS DE LA BASE DE DATOS*

% Ruta de la carpeta donde están las subcarpetas 'Grupo A' y 'Grupo B'
outputFolder = 'C:\Users\lolas\OneDrive\Documentos\Universidad\5º\TFG_Biomédica\BBDD\Grupos A y B_Punto\Grupos A y B_Punto\Grupos A y B';

% Verificar si la ruta existe
if exist(outputFolder, 'dir') == 0
    error('La ruta no existe. Verifica la ruta de la carpeta.');
end

% Verificar si hay archivos CSV en 'Grupo A'
grupoA_CSV = dir(fullfile(outputFolder, 'Grupo A', '*.csv'));
if isempty(grupoA_CSV)
    disp('No se encontraron archivos CSV en la carpeta Grupo A.');
else
    disp('Archivos CSV encontrados en Grupo A:');
    disp({grupoA_CSV.name});
end
% Verificar si hay archivos CSV en 'Grupo B'
grupoB_CSV = dir(fullfile(outputFolder, 'Grupo B', '*.csv'));
if isempty(grupoB_CSV)
    disp('No se encontraron archivos CSV en la carpeta Grupo B.');
else
    disp('Archivos CSV encontrados en Grupo B:');
    disp({grupoB_CSV.name});
end
%% 
% *GENERACIÓN DE DOS TABLAS PARA  ALMACENAR DATOS DE GRUPO A Y GRUPO B*
% 
% *GENERACIÓN TABLA CON DATOS DEL GRUPO A (PAIN-FREE)*

% Inicializar celdas para almacenar los datos de cada grupo
pacienteA_data = {};
pacienteB_data = {};

%DIVIDIR ENTRE 50 EL TIEMPO CON TICS

% Cargar archivos CSV de 'grupo A'
    csvPathA = fullfile(grupoA_CSV(1).folder, grupoA_CSV(1).name);
    opts = detectImportOptions(csvPathA, 'Delimiter', ';');
    opts.VariableNames = {'Sample', 'Time', 'Device', 'X', 'Y', 'Z', 'qX', 'qY', 'qZ', 'qW', 'x', 'y', 'z'};
    opts.VariableTypes = {'double', 'double', 'char', 'double', 'double', 'double', ...
                          'double', 'double', 'double', 'double', 'double', 'double', 'double'};
    pacienteA_data{1} = readtable(csvPathA,opts);
    fprintf('Paciente A - Archivo cargado: %s\n', grupoA_CSV(1).name);
% Paso a estructura de tabla
tablaPacienteA = pacienteA_data{1};
% Eliminar registros de calentamiento
idxStored = find(strcmp(tablaPacienteA.Device, 'stored.Apple 3'), 1, 'first');
if ~isempty(idxStored)
    tablaPacienteA = tablaPacienteA(idxStored:end, :);
else
     warning('No se encontró "stored.Apple 3" en el paciente %d', 1);
end
% Filtra los registros donde Device = 'H'
datosFiltradosA = tablaPacienteA(strcmp(tablaPacienteA.Device, 'H'), :);
disp(head(datosFiltradosA));
%% 
% *GENERACIÓN DE TABLA B CON DATOS DEL GRUPO B (PAIN-AFFECTED)*

% Cargar archivos CSV de 'grupo B'
    csvPathB = fullfile(grupoB_CSV(1).folder, grupoB_CSV(1).name);
    opts = detectImportOptions(csvPathB, 'Delimiter', ';');
    opts.VariableNames = {'Sample', 'Time', 'Device', 'X', 'Y', 'Z', 'qX', 'qY', 'qZ', 'qW', 'x', 'y', 'z'};
    opts.VariableTypes = {'double', 'double', 'char', 'double', 'double', 'double', ...
                          'double', 'double', 'double', 'double', 'double', 'double', 'double'};
    %opts.DataLines = [4, Inf];
    pacienteB_data{1} = readtable(csvPathB,opts);
    fprintf('Paciente B - Archivo cargado: %s\n', grupoB_CSV(1).name);
% Paso a estructura de tabla
tablaPacienteB = pacienteB_data{1};
% Eliminar registros de calentamiento
idxStored = find(strcmp(tablaPacienteB.Device, 'stored.Apple 3'), 1, 'first');
if ~isempty(idxStored)
    tablaPacienteB = tablaPacienteB(idxStored:end, :);
else
     warning('No se encontró "stored.Apple 3" en el paciente %d', 1);
end
% Filtra los registros donde Device = 'H'
datosFiltradosB = tablaPacienteB(strcmp(tablaPacienteB.Device, 'H'), :);
disp(head(datosFiltradosB));
%% 
% *ANÁLISIS DESPLAZAMIENTO*

% Recorremos cada paciente en el grupoA_data
    % Extraemos las columnas
    tiempo_A = datosFiltradosA(:, 2);
    posicionX_A = datosFiltradosA(:, 4);
    posicionY_A = datosFiltradosA(:, 5);
    posicionZ_A = datosFiltradosA(:, 6);
    anguloEulerX_A = datosFiltradosA(:, 11);
    anguloEulerY_A = datosFiltradosA(:, 12);
    anguloEulerZ_A = datosFiltradosA(:, 13);
    

    tiempo_B = datosFiltradosB(:, 2);
    posicionX_B = datosFiltradosB(:, 4);
    posicionY_B = datosFiltradosB(:, 5);
    posicionZ_B = datosFiltradosB(:, 6);
    anguloEulerX_B = datosFiltradosB(:, 11);
    anguloEulerY_B = datosFiltradosB(:, 12);
    anguloEulerZ_B = datosFiltradosB(:, 13);

% Visualizamos el resultado
%disp(posicionX);
%disp(posicionY);
%disp(posicionZ);
%% 
% *GRÁFICA DE LA POSICIÓN EN EJES X-Y-Z RESPECTO AL TIEMPO*

posicionX_A = table2array(posicionX_A);
posicionY_A = table2array(posicionY_A);
posicionZ_A = table2array(posicionZ_A);
tiempo_A = table2array(tiempo_A);
posicionX_B = table2array(posicionX_B);
posicionY_B = table2array(posicionY_B);
posicionZ_B = table2array(posicionZ_B);
tiempo_B = table2array(tiempo_B);
%% 
% 

%Centramos en 0 los datos--> restamos el promedio de las primeras 50 muestras
promedio_X_A = mean(posicionX_A(1:50));
promedio_Y_A = mean(posicionY_A(1:50));
promedio_Z_A = mean(posicionZ_A(1:50));

% Centrar las posiciones restando el promedio correspondiente
posicionX_A_centrada = posicionX_A - promedio_X_A;
posicionY_A_centrada = posicionY_A - promedio_Y_A;
posicionZ_A_centrada = posicionZ_A - promedio_Z_A;

%Centramos en 0 los datos--> restamos el promedio de las primeras 50 muestras
promedio_X_B = mean(posicionX_B(1:50));
promedio_Y_B = mean(posicionY_B(1:50));
promedio_Z_B = mean(posicionZ_B(1:50));

% Centrar las posiciones restando el promedio correspondiente
posicionX_B_centrada = posicionX_B - promedio_X_B;
posicionY_B_centrada = posicionY_B - promedio_Y_B;
posicionZ_B_centrada = posicionZ_B - promedio_Z_B;
%% 
% *GRÁFICAS GRUPO A POSICIÓN CENTRADAS*

figure;
% Gráfico de la posición X_A
subplot(3, 1, 1); % Primera subgráfica
plot(tiempo_A, posicionX_A_centrada, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('Position X');
title('Position X - Time (A)');
% Gráfico de la posición Y_A
subplot(3, 1, 2); % Segunda subgráfica
plot(tiempo_A, posicionY_A_centrada, 'g-', 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('Position Y');
title('Position Y - Time (A)');

% Gráfico de la posición Z_A
subplot(3, 1, 3); % Tercera subgráfica
plot(tiempo_A, posicionZ_A_centrada, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('Position Z');
title('Position Z - Time (A)');
%% 
% *GRÁFICAS GRUPO B POSICIÓN CENTRADAS*

figure;
% Gráfico de la posición X_B
subplot(3, 1, 1); % Primera subgráfica
plot(tiempo_B, posicionX_B_centrada, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('Position X');
title('Position X - Time (B)');
% Gráfico de la posición Y_B
subplot(3, 1, 2); % Segunda subgráfica
plot(tiempo_B, posicionY_B_centrada, 'g-', 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('Position Y');
title('Position Y - Time (B)');

% Gráfico de la posición Z_B
subplot(3, 1, 3); % Tercera subgráfica
plot(tiempo_B, posicionZ_B_centrada, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('Position Z');
title('Position Z - Time (B)');
%% 
% *ÁNGULOS EULER*

anguloEulerX_A = table2array(anguloEulerX_A);
anguloEulerY_A = table2array(anguloEulerY_A);
anguloEulerZ_A = table2array(anguloEulerZ_A);

%Centramos en 0 los datos--> restamos el promedio de las primeras 50 muestras
promedio_x_A = mean(anguloEulerX_A(1:50));
promedio_y_A = mean(anguloEulerY_A(1:50));
promedio_z_A = mean(anguloEulerZ_A(1:50));

% Centrar los ángulos restando el promedio correspondiente
anguloEulerX_A_centrado = anguloEulerX_A - promedio_x_A;
anguloEulerY_A_centrado = anguloEulerY_A - promedio_y_A;
anguloEulerZ_A_centrado = anguloEulerZ_A - promedio_z_A;

%Hago wrapping para que los ángulos que están entre [0,360] estén entre [-180,180] grados
anguloEulerX_A_centrado = wrapTo180(anguloEulerX_A_centrado);
anguloEulerY_A_centrado = wrapTo180(anguloEulerY_A_centrado);
anguloEulerZ_A_centrado = wrapTo180(anguloEulerZ_A_centrado);
%% 
% 
anguloEulerX_B = table2array(anguloEulerX_B);
anguloEulerY_B = table2array(anguloEulerY_B);
anguloEulerZ_B = table2array(anguloEulerZ_B);

%Centramos en 0 los datos--> restamos el promedio de las primeras 50 muestras
promedio_x_B = mean(anguloEulerX_B(1:50));
promedio_y_B = mean(anguloEulerY_B(1:50));
promedio_z_B = mean(anguloEulerZ_B(1:50));

% Centrar los ángulos restando el promedio correspondiente
anguloEulerX_B_centrado = anguloEulerX_B - promedio_x_B;
anguloEulerY_B_centrado = anguloEulerY_B - promedio_y_B;
anguloEulerZ_B_centrado = anguloEulerZ_B - promedio_z_B;

%Hago wrapping para que los ángulos que están entre [0,360] estén entre [-180,180] grados
anguloEulerX_B_centrado = wrapTo180(anguloEulerX_B_centrado);
anguloEulerY_B_centrado = wrapTo180(anguloEulerY_B_centrado);
anguloEulerZ_B_centrado = wrapTo180(anguloEulerZ_B_centrado);
%% 
% *GRÁFICA GRUPO A ÁNGULOS DE EULER CENTRADA* 

figure;
% Gráfico de x
subplot(3, 1, 1); % Primera subgráfica
plot(tiempo_A, anguloEulerX_A_centrado, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('Angle X');
title('Angle X - Time (A)');
% Gráfico de y
subplot(3, 1, 2); % Segunda subgráfica
plot(tiempo_A, anguloEulerY_A_centrado, 'g-', 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('Angle Y');
title('Angle Y - Time (A)');

% Gráfico de z
subplot(3, 1, 3); % Tercera subgráfica
plot(tiempo_A, anguloEulerZ_A_centrado, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('Angle Z');
title('Angle Z - Time (A)');
%% 
% *GRÁFICA GRUPO B ÁNGULOS DE EULER CENTRADA*

figure;
% Gráfico de x
subplot(3, 1, 1); % Primera subgráfica
plot(tiempo_B, anguloEulerX_B_centrado, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('Angle X');
title('Angle X - Time (B)');
% Gráfico de la posición Y
subplot(3, 1, 2); % Segunda subgráfica
plot(tiempo_B, anguloEulerY_B_centrado, 'g-', 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('Angle Y');
title('Angle Y - Time (B)');

% Gráfico de la posición Z
subplot(3, 1, 3); % Tercera subgráfica
plot(tiempo_B, anguloEulerZ_B_centrado, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('Angle Z');
title('Angle Z - Time (B)');
%% 
% *MÁXIMOS, MEDIA Y RANGO DE MOVIMIENTO EN CADA EJE (CENTRADO)*

fprintf('Máximo ángulo en el eje X del grupo A: %.2f grados\n',max(anguloEulerX_A_centrado));
fprintf('Máximo ángulo en el eje Y del grupo A: %.2f grados\n',max(anguloEulerY_A_centrado));
fprintf('Máximo ángulo en el eje Z del grupo A: %.2f grados\n',max(anguloEulerZ_A_centrado));
fprintf('Ángulo medio en el eje X del grupo A: %.2f grados\n',mean(anguloEulerX_A_centrado));
fprintf('Ángulo medio en el eje Y del grupo A: %.2f grados\n',mean(anguloEulerY_A_centrado));
fprintf('Ángulo medio en el eje Z del grupo A: %.2f grados\n',mean(anguloEulerZ_A_centrado));
ROM_X_A = max(anguloEulerX_A_centrado) - min(anguloEulerX_A_centrado);
ROM_Y_A = max(anguloEulerY_A_centrado) - min(anguloEulerY_A_centrado);
ROM_Z_A = max(anguloEulerZ_A_centrado) - min(anguloEulerZ_A_centrado); 
fprintf('Rango de movimiento en el eje X del grupo A (flexión/extensión): %.2f grados\n', ROM_X_A);
fprintf('Rango de movimiento en el eje Y del grupo A (inclinación lateral): %.2f grados\n', ROM_Y_A);
fprintf('Rango de movimiento en el eje Z del grupo A (rotación): %.2f grados\n', ROM_Z_A);
fprintf('Máximo ángulo en el eje X del grupo B: %.2f grados\n',max(anguloEulerX_B_centrado));
fprintf('Máximo ángulo en el eje Y del grupo B: %.2f grados\n',max(anguloEulerY_B_centrado));
fprintf('Máximo ángulo en el eje Z del grupo B: %.2f grados\n',max(anguloEulerZ_B_centrado));
fprintf('Ángulo medio en el eje X del grupo B: %.2f grados\n',mean(anguloEulerX_B_centrado));
fprintf('Ángulo medio en el eje Y del grupo B: %.2f grados\n',mean(anguloEulerY_B_centrado));
fprintf('Ángulo medio en el eje Z del grupo B: %.2f grados\n',mean(anguloEulerZ_B_centrado));
ROM_X_B = max(anguloEulerX_B_centrado) - min(anguloEulerX_B_centrado); 
ROM_Y_B = max(anguloEulerY_B_centrado) - min(anguloEulerY_B_centrado); 
ROM_Z_B = max(anguloEulerZ_B_centrado) - min(anguloEulerZ_B_centrado); 
fprintf('Rango de movimiento en el eje X del grupo B (flexión/extensión): %.2f grados\n', ROM_X_B);
fprintf('Rango de movimiento en el eje Y del grupo B (inclinación lateral): %.2f grados\n', ROM_Y_B);
fprintf('Rango de movimiento en el eje Z del grupo B (rotación): %.2f grados\n', ROM_Z_B);
%% 
% *CÁLCULO DE LAS VELOCIDADES (ver si hace falta filtrado)*

dt = (1/50); %Periodo de muestreo = 20ms
%La velocidad es la primera derivada
velocidad_X_A = diff(posicionX_A_centrada)/dt;
velocidad_Y_A = diff(posicionY_A_centrada)/dt;
velocidad_Z_A = diff(posicionZ_A_centrada)/dt;

velocidad_X_B = diff(posicionX_B_centrada)/dt;
velocidad_Y_B = diff(posicionY_B_centrada)/dt;
velocidad_Z_B = diff(posicionZ_B_centrada)/dt;

%Filtrado de las velocidades y aceleraciones
fs = 50;  % Frecuencia de muestreo (Hz)
fc = 20;  % Frecuencia de corte (Hz) 
orden = 2;  % Orden del filtro 
[b, a] = butter(orden, fc/(fs/2), 'low');

velocidad_X_A = filtfilt(b, a, velocidad_X_A);
velocidad_Y_A = filtfilt(b, a, velocidad_Y_A);
velocidad_Z_A = filtfilt(b, a, velocidad_Z_A);
aceleracion_X_A = diff(velocidad_X_A)/dt;
aceleracion_Y_A = diff(velocidad_Y_A)/dt;
aceleracion_Z_A = diff(velocidad_Z_A)/dt;

velocidad_X_B = filtfilt(b, a, velocidad_X_B);
velocidad_Y_B = filtfilt(b, a, velocidad_Y_B);
velocidad_Z_B = filtfilt(b, a, velocidad_Z_B);
aceleracion_X_B = diff(velocidad_X_B)/dt;
aceleracion_Y_B = diff(velocidad_Y_B)/dt;
aceleracion_Z_B = diff(velocidad_Z_B)/dt;
%% 
% *GRÁFICA DE LA VELOCIDAD RESPECTO AL TIEMPO PARA EL GRUPO A*

figure;
% Gráfico en X
subplot(3, 1, 1); 
plot(tiempo_A(1:end-1), velocidad_X_A, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('Speed X');
title('Speed X - Time (A)');
% Gráfico en Y
subplot(3, 1, 2); % Segunda subgráfica
plot(tiempo_A(1:end-1), velocidad_Y_A, 'g-', 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('Speed Y');
title('Speed Y - Time (A)');

% Gráfico en Z
subplot(3, 1, 3); % Tercera subgráfica
plot(tiempo_A(1:end-1), velocidad_Z_A, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('Speed Z');
title('Speed Z - Time (A)');
%% 
% *GRÁFICA DE LA VELOCIDAD RESPECTO AL TIEMPO PARA EL GRUPO B*

figure;
% Gráfico en X
subplot(3, 1, 1); 
plot(tiempo_B(1:end-1), velocidad_X_B, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('Speed X');
title('Speed X - Time (B)');
% Gráfico en Y
subplot(3, 1, 2); % Segunda subgráfica
plot(tiempo_B(1:end-1), velocidad_Y_B, 'g-', 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('Speed Y');
title('Speed Y - Time (B)');

% Gráfico en Z
subplot(3, 1, 3); % Tercera subgráfica
plot(tiempo_B(1:end-1), velocidad_Z_B, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('Speed Z');
title('Speed Z - Time (B)');
%% 
% *GRÁFICA DE LA ACELERACIÓN RESPECTO AL TIEMPO PARA EL GRUPO A*

figure;
% Gráfico en X
subplot(3, 1, 1); 
plot(tiempo_A(1:end-2), aceleracion_X_A, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('Acceleration X');
title('Acceleration X - Time (A)');
% Gráfico en Y
subplot(3, 1, 2); % Segunda subgráfica
plot(tiempo_A(1:end-2), aceleracion_Y_A, 'g-', 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('Acceleration Y');
title('Acceleration Y - Time (A)');

% Gráfico en Z
subplot(3, 1, 3); % Tercera subgráfica
plot(tiempo_A(1:end-2), aceleracion_Z_A, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('Acceleration Z');
title('Acceleration Z - Time (A)');
%% 
% *GRÁFICA DE LA ACELERACIÓN RESPECTO AL TIEMPO PARA EL GRUPO B*

figure;
% Gráfico en X
subplot(3, 1, 1); 
plot(tiempo_B(1:end-2), aceleracion_X_B, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('Acceleration X');
title('Acceleration X - Time (B)');
% Gráfico en Y
subplot(3, 1, 2); % Segunda subgráfica
plot(tiempo_B(1:end-2), aceleracion_Y_B, 'g-', 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('Acceleration Y');
title('Acceleration Y - Time (B)');

% Gráfico en Z
subplot(3, 1, 3); % Tercera subgráfica
plot(tiempo_B(1:end-2), aceleracion_Z_B, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('Acceleration Z');
title('Acceleration Z - Time (B)');
