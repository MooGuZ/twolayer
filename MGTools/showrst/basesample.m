% This script generate handmade bases to show the discriptive power of
% complex bases

sz    = 256;
[x,y] = meshgrid(1:sz,1:sz);
z     = complex(x,y); clear x y

ncycle = 7;

%% Purely Rotation
center = complex((1+sz)/2,(1+sz)/2);
phase = angle(z - center);
base = exp(1j*phase(:)*ncycle);

%% Water Flow
c1 = complex(1,64); c2 = complex(sz,sz-64);
base = complex(zeros(sz,sz));
dist = abs([z(:)-c1,z(:)-c2]);
[~,index] = min(dist,[],2);
base(index==1) = z(index==1) - c1;
base(index==2) = z(index==2) - c2;
base = angle(base);
base(index==1) = angle(c2-c1) - base(index==1);
base(index==2) = base(index==2) - angle(c1-c2);
base = exp(1j*base(:)*ncycle);

%% Refraction
c1 = complex(64,1); c2 = complex(sz-64,sz);
phase = zeros(sz,sz); amplitude = zeros(sz,sz);
index = angle(z-c1) < angle(c2-c1);
phase(index) = angle(z(index)-c1) - angle(c2-c1);
phase(~index) = angle(c1-c2) - angle(z(~index)-c2);
dist = abs([z(:)-c1,z(:)-c2]); d = abs(c1-c2);
amplitude(index(:) & (dist(:,1) > d/3) & (dist(:,1) < 2*d/3)) = 1;
amplitude(~index(:) & (dist(:,2) > d/3) & (dist(:,2) < 2*d/3)) = 1;
H = fspecial('gaussian',7,3);
amplitude = imfilter(amplitude,H,'replicate');
base = amplitude(:) .* exp(1j*phase(:)*ncycle);

%% Scaling
center = complex((1+sz)/2,(1+sz)/2);
phase = 2*pi*abs(z-center)/sz;
base = exp(1j*phase(:)*ncycle);

%% Shifting
phase = 2*pi*real(z)/sz;
base = exp(1j*phase(:)*ncycle);