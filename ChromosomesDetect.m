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

%%
I_green = imread("Chromosomes-green.tif") ;
% figure ; imshow(I_green) ;
% imhist(I_green)

% Method 1 : Otsu 

[T,EM] = graythresh(I_green);
Iotsu = I_green; 
Iotsu(Iotsu>T*(2^16 -1)) = 2^16 -1;
Iotsu(Iotsu<=T*(2^16 -1)) = 0;
Iotsu_green = logical(Iotsu);

%figure ; imshow(logical(Iotsu_green));

% Method 2 : manual threshold
I_green_threshold = I_green;
I_green_threshold(I_green_threshold > 15000) = 2^16 ;
I_green_threshold(I_green_threshold <= 15000) = 0 ;

% figure ; imshow(logical(I_green_threshold))

cells = cell(zones, 1) ;
coord = zeros(2,2,zones) ;

for z=1:zones 
    [r, c] = find(L==z);
    burden = zeros(max(r)-min(r)+1, max(c)-min(c)+1);
    for i=min(r):max(r)
        for j=min(c):max(c)
            if L(i,j) == z && I_blue_threshold1(i,j)==2^16-1
                burden(i-min(r)+1,j-min(c)+1) = 1;
            else 
                burden(i-min(r)+1,j-min(c)+1) = 0;
            end
        end
    end
    burden = logical(burden);
    [r, c] = find(burden==1);
    coord(:, :, z) = [min(r) max(r) ; min(c) max(c)];
    burden = burden(min(r):max(r), min(c):max(c)) ;
    cells{z} = burden ;
end

% Plot the cells
figure;
grid = double(4);
q = double(mod(zones, grid));
if q > 0
    p = double(zones-q);
    for plotId=1:p
        subplot(grid, p/grid, plotId) ;
        imshow(cells{plotId}) ;
    end
    figure;
    for plotId=1:q
        subplot(1, q, plotId) ;
        imshow(cells{zones-plotId+1}) ;
    end
else
    for plotId=1:zones
        subplot(grid, grid, plotId) ;
        imshow(cells{plotId}) ;
    end
end


