%% Test 3
 
%%Alternative
[a1, aFs1] = audioread('AL1.mp3');
[a2, aFs2] = audioread('AL2.mp3');
[a3, aFs3] = audioread('AL3.m4a');
[a4, aFs4] = audioread('AL4.m4a');
[a5, aFs5] = audioread('AL5.mp3');
 
%%Classical
[b1, bFs1] = audioread('CL1.mp3');
[b2, bFs2] = audioread('CL2.mp3');
[b3, bFs3] = audioread('CL3.mp3');
[b4, bFs4] = audioread('CL4.mp3');
[b5, bFs5] = audioread('CL5.mp3');
 
%%Hip-Hop
[c1, cFs1] = audioread('HH1.mp3');
[c2, cFs2] = audioread('HH2.m4a');
[c3, cFs3] = audioread('HH3.mp3');
[c4, cFs4] = audioread('HH4.mp3');
[c5, cFs5] = audioread('HH5.mp3');
%%
Fs = aFs1;
 
a1com = (a1(:, 1) + a1(:, 2))/2;
a1Seg = a1com(Fs*60 + 1:Fs*65);
a2com = (a2(:, 1) + a2(:, 2))/2;
a2Seg = a2com(Fs*60 + 1:Fs*65);
a3com = (a3(:, 1) + a3(:, 2))/2;
a3Seg = a3com(Fs*60 + 1:Fs*65);
a4com = (a4(:, 1) + a4(:, 2))/2;
a4Seg = a4com(Fs*60 + 1:Fs*65);
a5com = (a5(:, 1) + a5(:, 2))/2;
a5Seg = a5com(Fs*60 + 1:Fs*65);
a6Seg = a1com(Fs*66 + 1:Fs*71);
a7Seg = a2com(Fs*66 + 1:Fs*71);
a8Seg = a3com(Fs*66 + 1:Fs*71);
a9Seg = a4com(Fs*66 + 1:Fs*71);
a10Seg = a5com(Fs*66 + 1:Fs*71);
 
b1com = (b1(:, 1) + b1(:, 2))/2;
b1Seg = b1com(Fs*60 + 1:Fs*65);
b2com = (b2(:, 1) + b2(:, 2))/2;
b2Seg = b2com(Fs*60 + 1:Fs*65);
b3com = (b3(:, 1) + b3(:, 2))/2;
b3Seg = b3com(Fs*60 + 1:Fs*65);
b4com = (b4(:, 1) + b4(:, 2))/2;
b4Seg = b4com(Fs*60 + 1:Fs*65);
b5com = (b5(:, 1) + b5(:, 2))/2;
b5Seg = b5com(Fs*60 + 1:Fs*65);
b6Seg = b1com(Fs*66 + 1:Fs*71);
b7Seg = b2com(Fs*66 + 1:Fs*71);
b8Seg = b3com(Fs*66 + 1:Fs*71);
b9Seg = b4com(Fs*66 + 1:Fs*71);
b10Seg = b5com(Fs*66 + 1:Fs*71);
 
 
c1com = (c1(:, 1) + c1(:, 2))/2;
c1Seg = c1com(Fs*60 + 1:Fs*65);
c2com = (c2(:, 1) + c2(:, 2))/2;
c2Seg = c2com(Fs*60 + 1:Fs*65);
c3com = (c3(:, 1) + c3(:, 2))/2;
c3Seg = c3com(Fs*60 + 1:Fs*65);
c4com = (c4(:, 1) + c4(:, 2))/2;
c4Seg = c4com(Fs*60 + 1:Fs*65);
c5com = (c5(:, 1) + c5(:, 2))/2;
c5Seg = c5com(Fs*60 + 1:Fs*65);
c6Seg = c1com(Fs*66 + 1:Fs*71);
c7Seg = c2com(Fs*66 + 1:Fs*71);
c8Seg = c3com(Fs*66 + 1:Fs*71);
c9Seg = c4com(Fs*66 + 1:Fs*71);
c10Seg = c5com(Fs*66 + 1:Fs*71);
%%
aSeg = [a1Seg' ; a2Seg'; a3Seg'; a4Seg'; a5Seg'; a6Seg'; a7Seg'; a8Seg'; a9Seg'; a10Seg'];
bSeg = [b1Seg' ; b2Seg'; b3Seg'; b4Seg'; b5Seg'; b6Seg'; b7Seg'; b8Seg'; b9Seg'; b10Seg'];
cSeg = [c1Seg' ; c2Seg'; c3Seg'; c4Seg'; c5Seg'; c6Seg'; c7Seg'; c8Seg'; c9Seg'; c10Seg'];
 
%%
L = 10;
n = 5*Fs;
rate = 0.1;
tslide = 0:rate:5;
window = -20;
%% Windowed Fourier Transform Plot
 ai = aSeg(1, :);
 aSpec = [];
 for j = 1:length(tslide)
    	t = (1:length(ai))/Fs;
    	g = exp(window*(t-tslide(j)).^2); %Guassian Filter
    	aG = g.*ai;
    	aG_fft = fft(aG); 
    	aSpec = [aSpec; resample(abs(fftshift(aG_fft)), 1, 10)];
    	subplot(3, 1, 1);
    	plot(t, ai,'k',t, g,'r');
    	subplot(3, 1, 2);
    	plot(t, aG);
    	n = 220500;
    	k = (2*pi)/5*[0:n/2-1 -n/2:-1];
    	ks = fftshift(k);
    	subplot(3, 1, 3);
    	aGt = aG_fft(1, n/2:n);
    	plot(ks(n/2:n), abs(fftshift(aGt)/max(abs(aGt))));
    	drawnow;
	end
%%
xa = [];
for i = 1:L
	ai = aSeg(i, :);
	aSpec = [];
	for j = 1:length(tslide)
    	t = (1:length(ai))/Fs;
    	g = exp(window*(t-tslide(j)).^2); %Guassian Filter
    	aG = g.*ai;
    	aG_fft = fft(aG); 
    	aSpec = [aSpec; resample(abs(fftshift(aG_fft)), 1, 10)];
	end
	xa = [xa, aSpec];
end
%% pca for a
[m, n] = size(xa);
mn = mean(xa, 2);
xa_sub=xa-repmat(mn, 1, n);
[u, s, v] = svd(xa_sub/sqrt(n-1), 'econ');
y = u'*xa_sub;
Cya = y*y.'/(n-1);
PCa = diag(Cya);
save PCa.mat PCa;
figure (1);
subplot(1, 3, 1);
scatter(1:length(PCa), PCa);
title('PCs for Alternative');
xlabel('Component #');
ylabel('Singular Value');
%%
xb = [];
for i = 1:L
	bi = bSeg(i, :);
	bSpec = [];
	for j = 1:length(tslide)
    	t = (1:length(bi))/Fs;
    	g = exp(window*(t-tslide(j)).^2); %Guassian Filter
    	bG = g.*bi;
    	bG_fft = fft(bG); 
    	bSpec = [bSpec; resample(abs(fftshift(bG_fft)), 1, 10)];
	end
	xb = [xb, bSpec];
end
%% pca for b
[m, n] = size(xb);
mn = mean(xb, 2);
xb_sub=xb-repmat(mn, 1, n);
[u, s, v] = svd(xb_sub/sqrt(n-1), 'econ');
y = u'*xb_sub;
Cyb = y*y.'/(n-1);
PCb = diag(Cyb);
save PCb.mat PCb;
subplot(1, 3, 2);
scatter(1:length(PCb), PCb);
title('Classical');
xlabel('Component #');
ylabel('Singular Value');
%%
xc = [];
for i = 1:L
	ci = cSeg(i, :);
	cSpec = [];
	for j = 1:length(tslide)
    	t = (1:length(ci))/Fs;
    	g = exp(window*(t-tslide(j)).^2); %Guassian Filter
    	cG = g.*ci;
    	cG_fft = fft(cG); 
    	cSpec = [cSpec; resample(abs(fftshift(cG_fft)), 1, 10)];
	end
	xc = [xc, cSpec];
end
%% pca for c
[m, n] = size(xc);
mn = mean(xc, 2);
xc_sub=xc-repmat(mn, 1, n);
[u, s, v] = svd(xc_sub/sqrt(n-1), 'econ');
y = u'*xc_sub;
Cyc = y*y.'/(n-1);
PCc = diag(Cyc);
save PCc.mat PCc;
subplot(1, 3, 3);
scatter(1:length(PCc), PCc);
title('Hip-Hop');
xlabel('Component #');
ylabel('Singular Value');
%% pca
x = [xa, xb, xc];
[m, n] = size(x);
mn = mean(x, 2);
x_sub=x-repmat(mn, 1, n);
[u, s, v] = svd(x_sub/sqrt(n-1), 'econ');
y = u'*x_sub;
 
%%
rand = randperm (L);
x1 = y(1:14, 1:L);
x2 = y(1:14, L+1:2*L);
x3 = y(1:14, 2*L+1:3*L);
 
xtrain = [x1(:, rand(1:6)), x2(:, rand(1:6)), x3(:, rand(1:6))];
xtest = [x1(:, rand(7:10)), x2(:,rand(7:10)), x3(:, rand(7:10))];
trainInto = [ones(1, 6), 2*ones(1, 6), 3*ones(1, 6)];
 
result = classify(xtest', xtrain', trainInto);
figure(2);
bar(result);
title('Classification');
 
