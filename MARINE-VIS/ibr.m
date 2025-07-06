function img_ibr = ibr(originalImg, r_min_pct, r_max_pct)
    % Validación de parámetros
    if nargin < 2, r_min_pct = 0.0027; end  % 0.27% por defecto
    if nargin < 3, r_max_pct = 0.35; end     % 35% por defecto
    
    % 1. Conversión a espacio HSV
    hsv_img = rgb2hsv(originalImg);
    V_original = hsv_img(:,:,3);
    
    % 2. Calcular parámetros basados en el lado más corto
    [rows, cols] = size(V_original);
    shorter = min([rows, cols]);
    r_min = shorter * r_min_pct;       % Radio mínimo para iluminación
    r_max = shorter * r_max_pct;       % Radio máximo para iluminación

    % 3. Transformada de Fourier
    F = fft2(double(V_original));
    Fshift = fftshift(F);
    
    % 4. Crear sistema de coordenadas espectrales
    [X, Y] = meshgrid(1:cols, 1:rows);
    center_x = ceil(cols/2);
    center_y = ceil(rows/2);
    D = sqrt((X - center_x).^2 + (Y - center_y).^2);
    
    % 5. Crear máscara 
    % - Máscara pasa-banda (corrección de iluminación)
    mask_CI = (D >= r_min) & (D <= r_max);
    
    % 6. Aplicar máscara al espectro
    Fshift_filtered = Fshift .* mask_CI;
    
    % 7. Transformada inversa de Fourier
    V_corrected = real(ifft2(ifftshift(Fshift_filtered)));
    V_corrected = mat2gray(V_corrected);
    
    % 8. Reconstrucción de la imagen RGB
    hsv_img(:,:,3) = V_corrected;
    img_ibr = hsv2rgb(hsv_img);

    %9 Normalización del resultado
    img_ibr = mat2gray(img_ibr);
end