%----------------------------------------------------------------------
%                      Frequency Domain Decomposition 
%----------------------------------------------------------------------
%   
%   The algorithm 'fdd_k' identifies modal parameters in the frequency
%   domain using frequency domain decomposition (FDD) and peak picking.
%
%                 [Frq,phi]=fdd_k(y,fs,nfft,w)
% 
%   Inputs:
%           y:    matrix of measured outputs as column vectors
%           fs:   sample frequency
%           nfft: number of nfft points
%           W:    window function 
%                   for example:
%                      hann:      hanning window 
%                      hamming:   hamming window
%                      chebwin:   Chebyshev window
%                      rectwin:   rectangular window
%           
%   Outputs:
%            Frq: natural frequency
%            phi: mode shape
%
% This algorithm was developed from the following literature:
% 1. Brincker, R. Ventura, C. (2015) Introduction to operational Modal Analysis, Wiley.
% 2. Brandt, A. (2011) Noise and Vibration Analysis Signal Analysis and Experimental Procedures, John Wiley and Sons.

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                         Keith Soal 2016
%                 questions to keithsoal at gmail.com
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function [Frq,phi]=fdd_k(y,fs,nfft,w)

% *********************************************
%               Measurement parameters 
% *********************************************
delta_freq = fs/nfft;   % Frequency resolution
NFFT=nfft/2+1;
sy = size(y);

% **************************************************************
%               Power spectral density matrix 
% **************************************************************
Gxy = zeros(sy(2),sy(2),NFFT);
F = zeros(sy(2),sy(2),NFFT);
% efficient cross spectral density estimation
[Gxy,f] = cpsd_k(y,w,nfft,fs);

% **************************************************************
%               Singular value decomposition 
% **************************************************************
S_n = input('Number of SVD plots ');      
for k = 1:size(Gxy,3)
    [U,S,V] = svd(Gxy(:,:,k));
        for i = 1:S_n
            Sv(k,i) = S(i,i); % Singular values
        end
    v(:,k) = U(:,1); % Mode shapes
end

% Singular value plotting
figure()
color_def = [0 0 1;0 1 0;1 0 0];
magnta = [1 0 1];
for nn=1:S_n
    if nn <= 3
        hold on; plot(f,mag2db(Sv(:,nn)),'color',color_def(nn,:));
        title('Singular Values')
        xlabel('Frequency [Hz]')
        ylabel('s amplitude')
        grid on
    else
        hold on; plot(f,mag2db(Sv(:,nn)),'color',magnta);
    end
end 

% ***********************************************************************
%               Modal paramter identification 
% ***********************************************************************
Fp=[];% Frequencies related to selected peaks

% **********************************************
%               Peak picking 
% **********************************************
NumPeaks=input('Enter the number of desired peaks:');
display('-----------------------------------------------------')
display('Peak selection procedure')
display('a: Draw rectangles around peaks')
display('b: Press "Space" key to continue the peak selection')
display('c: Press "any other key" if you have selected a peak by mistake and want to ignore it')
k=0;
while k~=NumPeaks
    A=getrect;  % Draw a rectangle around the peak
    [~,P1]=min(abs(f-A(1)));
    [~,P2]=min(abs(f-(A(1)+A(3))));
    [~,B]=max(Sv(P1:P2,1));
    Max=B+P1-1; % Frequency at the selected peak
    scatter(f(Max),mag2db(Sv(Max,1)),'MarkerEdgeColor','b','MarkerFaceColor','b') % Mark this peak
    pause;key=get(gcf,'CurrentKey');
    Fp(end+1,:)=[Max,f(Max)];
    if strcmp(key,'space')
        % Press space to continue peak selection
        k=k+1;
        scatter(f(Max),mag2db(Sv(Max,1)),'MarkerEdgeColor','g','MarkerFaceColor','g') % Mark this peak as green
    else
        % Press any other key to ignore this peak
        Fp(end,:)=[];
        scatter(f(Max),mag2db(Sv(Max,1)),'MarkerEdgeColor','r','MarkerFaceColor','r') % Mark this peak as red
    end
end
% Number selected peaks, respectively
[~,Sr]=sort(Fp(:,2));
Fp=Fp(Sr,:);

% Identified modal frequencies
Frq=Fp(:,2);

% **********************************************
%               Mode shapes
% **********************************************
for i=1:size(Fp,1)
    [Um, ~, ~] = svd(Gxy(:,:,Fp(i,1)));
    phi(:,i) = Um(:,1);
end

% Eigen vector plotting
figure()
for nn=1:size(Fp,1)
    subplot(size(Fp,1),1,nn)
    hold on; plot(real(phi(:,nn)),'b');
    title(['Mode Shape' num2str(nn)])
    xlabel('node point')
    ylabel('amplitude')
    grid on
    axis tight
end 

% **********************************************
%               Complexity plot
% **********************************************
figure()
rp = real(phi);
ip = imag(phi);
for nn=1:size(Fp,1)
    subplot(ceil(size(Fp,1)/2),2,nn)
    for i = 1:sy(2)
    hold on; plotv([rp(i,nn);ip(i,nn)],'b');
    end
    title(['Complexity plot mode' num2str(nn)])
    grid on
    axis equal
end 

% **********************************************
%                     MAC
% **********************************************
Mac = zeros(size(Fp,1),size(Fp,1));
for i=1:size(phi,2)
    for j=1:size(phi,2)
        Mac(i,j)=mac(phi(:,i),phi(:,j));
    end
end

figure()
M = size(phi,2); % grid resolution
x = linspace(1,M,M); % x-grid
y = linspace(1,M,M); % y-grid
bar3c(Mac,x,y,[],jet(50))
colorbar
colormap(jet(50))

end



