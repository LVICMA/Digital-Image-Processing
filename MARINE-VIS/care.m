function img_care = care(img_lbr, scales, alpha, beta, Gain, offset, s1, s2)
    if nargin < 2, scales = [15, 80, 250]; end
    if nargin < 3, alpha = 125; end
    if nargin < 4, beta = 46; end
    if nargin < 5, Gain = 30; end
    if nargin < 6, offset = -6; end
    if nargin < 7, s1 = 1; end
    if nargin < 8, s2 = 1; end

    img_care = zeros(size(img_lbr));
    msrcr();
    
    %% Chrominance Adjustment and Restoration Enhancement (CARE)
    % Técnicas utilizadas:
    % * MSRCR de Petro y otros, 2014
    % * Simplest Color Balance
    % * Normalización Lineal
    % II. Equilibrar la crominancia de la imagen para mejorar su contraste.
    function msrcr()
        [M, N, C] = size(img_lbr);        
        msr = zeros(M, N);  % Inicializar acumulador
            
        % Paso 1: Multi-Scale Retinex (MSR)
        for si = 1:numel(scales)
            sigma = scales(si);
            % Aplicar filtro Gaussiano
            v = gaussian_kernel(img_lbr, sigma);
            % Aplicar SSR
            diff = log(img_lbr + 0.001) - log(v + 0.001);
            msr = msr + (1/numel(scales)) * diff;   
        end

        % Paso 2: Restauración de Color (MSRCR)
        CRF = beta * (log(alpha * img_lbr) - log(sum(img_lbr,3))); 
        msrcr = exp(msr) .* CRF;
         
        % Paso 3: Simplest Color Balance
        out = SimplestColorBalance(msrcr, s1, s2);
        
        % Paso 4: Post-procesamiento
        img_care =  Gain * (out - offset);
        img_care = (img_care - min(img_care(:))) ./ (max(img_care(:)) - min(img_care(:)));

        function v = gaussian_kernel(I, sigma)
            [M, N, C] = size(I);  % Dimensiones originales y canales
            v = zeros(M, N, C);    % Inicializar imagen de salida
    
            for i = 1:C  % Procesar cada canal RGB
                u = I(:,:,i);
        
                % 1. Simetrización (extender a 2M x 2N)
                U = [u fliplr(u); flipud(u) rot90(u, 2)];
        
                % 2. FFT de la imagen simetrizada
                U_hat = fft2(U);
        
                % 3. Generar kernel gaussiano en dominio de frecuencia
                [P, Q] = size(U);
                [Ky, Kx] = meshgrid(0:Q-1, 0:P-1);
                omega_x = (2 * pi * Kx) / P;  % Frecuencias en eje x
                omega_y = (2 * pi * Ky) / Q;  % Frecuencias en eje y
                gaussian_ft = exp(-sigma^2/2 * (omega_x.^2 + omega_y.^2));
        
                % 4. Multiplicación en frecuencia
                V_hat = U_hat .* gaussian_ft;
        
                % 5. IFFT y obtener parte real
                V = real(ifft2(V_hat));
        
                % 6. Recortar al tamaño original
                v(:,:,i) = V(1:M, 1:N);
            end
        end
        function outImg = SimplestColorBalance(img, s1, s2)
            % Guardar la clase original de la imagen
            originalClass = class(img);
    
            % Convertir a double si es entero (mantiene rango 0-255)
            img = double(real(img));

            [M, N, C] = size(img);
            outImg = zeros(M, N, C, 'like', img); % Mantener tipo temporal
    
            % Valores precalculados
            N_pixels = M * N;
            bins = 0:256; % Bins para histograma
    
            for j = 1:C
                canal = img(:,:,j);
        
                % Histograma y acumulado
                counts = histcounts(canal, bins);
                histo_cum = cumsum(counts);
        
                % Cálculo de vmin
                if s1 == 0
                    vmin = min(canal(:));
                else
                    target_low = N_pixels * s1 / 100;
                    idx = find(histo_cum >= target_low, 1);
                    if isempty(idx)
                        vmin = 255;
                    else
                        vmin = idx - 1;
                    end
                end
        
                % Cálculo de vmax con ajuste vectorizado
                if s2 == 0
                    vmax = max(canal(:));
                else
                    target_high = N_pixels * (1 - s2/100);
                    idx = find(histo_cum <= target_high, 1, 'last');
                    if isempty(idx)
                        vmax = 0;
                    else
                        vmax = idx - 1;
                    end
                    vmax = min(vmax + (vmax < 254), 255); % Ajuste vectorizado
                end
        
                % Escalado y saturación
                canal_sat = min(max(canal, vmin), vmax);
                if vmax ~= vmin
                    canal_sat = (canal_sat - vmin) * (255 / (vmax - vmin));
                end
        
                outImg(:,:,j) = canal_sat;
            end
    
            % Conversión final al tipo original
            outImg = cast(outImg, originalClass);
        end
    end
end