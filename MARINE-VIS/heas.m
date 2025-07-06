function img_heas = heas(img_care)
    img_heas = zeros(size(img_care));
    segmentadorhibrido();

    %% Segmentación híbrida (Sobel + Otsu por refinamiento iterativo)
    % III. Segmentar la tilapia del fondo
    function segmentadorhibrido()
        % 1. Conversión eficiente a escala de grises
        img_gris = 0.299 * img_care(:,:,1) + 0.587 * img_care(:,:,2) + 0.114 * img_care(:,:,3);
    
        % 2. Detección inicial de bordes
        bordes = edge(img_gris, 'sobel', 0.05);
    
        % 3. Umbralización Otsu inicial
        mascara_otsu = imbinarize(img_gris, graythresh(img_gris));
    
        % 4. Refinamiento iterativo
        for i = 1:3
            % 4.1 Combinar bordes con máscara actual
            mascara_otsu = mascara_otsu & ~bordes;
        
            % 4.2 Operaciones morfológicas adaptativas
            se = crearElementoEstructurante(img_gris);
            mascara_otsu = imclose(mascara_otsu, se);
        
            % 4.3 Actualizar detección de bordes
            bordes = edge(img_gris, 'sobel', 0.05 + 0.02*i);
        end
        
        % 5. Máscara segmentada
        img_heas = imfill(mascara_otsu, 'holes');
    end
    %% Función auxiliar para crear elemento estructurante orientado
    function se = crearElementoEstructurante(img_gris)
        [alto, ancho] = size(img_gris);
        relacion_aspecto = ancho/alto;
    
        if relacion_aspecto > 1 % Imagen más ancha que alta
            se = strel('rectangle', [1 round(ancho/50)]);
        else
            se = strel('rectangle', [round(alto/50) 1]);
        end
    end
end