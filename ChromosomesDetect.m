I_blue = imread("Chromosomes-blue.tif") ;

% figure ;
% subplot(1,2,1); imshow(I_green);

% Threshold : method 1 with histogram equalization + median filter + manual threshold:

I_blue_eq = histeq(I_blue);
I_blue_eq_med = medfilt2(I_blue_eq, [30 30]);

% figure;
% subplot(1,2,1); imshow(I_blue_eq);
% subplot(1,2,2); imshow(I_blue_eq_med);

% imhist(I_blue_eq_med);

I_blue_threshold1 = I_blue_eq_med ;
I_blue_threshold1(I_blue_threshold1 > 6.2*10^4) = 2^16 ;
I_blue_threshold1(I_blue_threshold1 <= 6.2*10^4) = 0 ;

% Threshold : method 2 only manual threshold :

% imhist(I_blue);
I_blue_threshold2 = I_blue;
I_blue_threshold2(I_blue_threshold2 > 5000) = 2^16 ;
I_blue_threshold2(I_blue_threshold2 <= 5000) = 0 ;
 

% Threshold : method 3 FFT :
F = fftshift(fft2(I_blue));
Fb = F ;

cote = 100 ;
x = size(F, 1) /2;
y = size(F, 2) /2;

for i = 1:size(F, 1)
    for j= 1:size(F, 2)
        if i > x+cote || i < x-cote || j > y+cote || j< y-cote
            Fb(i,j) = 0;
        end
    end
end

% Fb(x-c:x+c, y-c:y+c) = 0; high filter pass
rI = ifft2(ifftshift(Fb));

% figure ;
% subplot(2,2,1); imagesc(log(1+abs(F)));
% subplot(2,2,2); imagesc(log(1+abs(Fb)));
% subplot(2,2,3); imshow(I_blue);
% subplot(2,2,4); imshow(uint16(real(rI)));

% Threshold : method 4 Otsu :

[T,EM] = graythresh(I_blue_eq_med);
Iotsu = I_blue_eq_med; 
Iotsu(Iotsu>T*(2^16 -1)) = 2^16 -1;
Iotsu(Iotsu<=T*(2^16 -1)) = 0;
Iotsu_blue = logical(Iotsu);

% Comparisons :

% figure ;
% subplot(1,3,1); imshow(logical(I_blue_threshold1));
% subplot(1,3,2); imshow(logical(I_blue_threshold2));
% subplot(1,3,3); imshow(logical(Iotsu_blue));


% Watershed : detect cells postions 
map = bwdist(I_blue_threshold1);
L = watershed(map);
disp = labeloverlay(I_blue, L);

% figure ;
% subplot(1,2,1); imshow(logical(I_blue_threshold1));
% subplot(1,2,2); imshow(disp);

zones = max(unique(L)); % 13 distinct regions with threshold1

% Crop the cells from the image:

cells_blue = cell(zones, 1) ;
coord = zeros(2,2,zones) ;

for z=1:zones 
    display(z) 
    [r1, c1] = find(L==z);
    %coord(:, :, z) = [min(r1) max(r1) ; min(c1) max(c1)];
    burden = zeros(max(r1)-min(r1)+1, max(c1)-min(c1)+1);
    for i=min(r1):max(r1)
        for j=min(c1):max(c1)
            if L(i,j) == z && I_blue_threshold1(i,j)==2^16-1
                burden(i-min(r1)+1,j-min(c1)+1) = 1;
            else 
                burden(i-min(r1)+1,j-min(c1)+1) = 0;
            end
        end
    end
    burden = logical(burden);
    [r2, c2] = find(burden==1);
    % attention au changement de base / d'indice avec la matrice burden (+min(r1)-1 ; +min(c1)-1) :
    coord(:, :, z) = [min(r2)+min(r1)-1 max(r2)+min(r1)-1 ; min(c2)+min(c1)-1 max(c2)+min(c1)-1];
    burden = burden(min(r2):max(r2), min(c2):max(c2)) ;
    cells_blue{z} = burden ;
    clear r1 r2 c1 c2 ;
end

% Plot the cells from blue image
grid = double(4);
q = double(mod(zones, grid));
if q > 0
    figure;
    p = double(zones-q);
    for plotId=1:p
        subplot(grid, p/grid, plotId) ;
        imshow(cells_blue{plotId}) ;
    end
    figure;
    for plotId=1:q
        subplot(1, q, plotId) ;
        imshow(cells_blue{zones-plotId+1}) ;
    end
else
    figure;
    for plotId=1:zones
        subplot(grid, grid, plotId) ;
        imshow(cells_blue{plotId}) ;
    end
end

%% 
I_green = imread("Chromosomes-green.tif") ;
% figure ; imshow(I_green) ;
% imhist(I_green)

cells_green = cell(zones, 1);
for z=1:size(coord,3)
    xmin = coord(1,1,z);
    xmax = coord(1,2,z);
    ymin = coord(2,1,z);
    ymax = coord(2,2,z);
    cells_green{z} = I_blue(xmin:xmax, ymin:ymax);
    clear xmin xmax ymin ymax
end

% Plot the cells from blue image
grid = double(4);
q = double(mod(zones, grid));
if q > 0
    figure;
    p = double(zones-q);
    for plotId=1:p
        subplot(grid, p/grid, plotId) ;
        imshow(cells_green{plotId}) ;
    end
    figure;
    for plotId=1:q
        subplot(1, q, plotId) ;
        imshow(cells_green{zones-plotId+1}) ;
    end
else
    figure;
    for plotId=1:zones
        subplot(grid, grid, plotId) ;
        imshow(cells_green{plotId}) ;
    end
end
