function generateWaveSurface()
% GENERATEWAVESURFACE  Create a time-varying Gaussian surface and save results.
%
%   Generates a sinusoidally oscillating Gaussian surface
%   over a fixed mesh grid. The amplitude varies with time, producing a
%   "breathing" wave. Results are saved as:
%       - HDF5 file (wave.h5)
%       - GIF animation (wave.gif)
%       - AVI video (wave.avi)
%
%   Required files:
%       - surfaceData1200.mat   (optional, for context)
%       - surfMesh.mat          (must contain xMesh and yMesh)
%
%   Output files:
%       - wave.h5
%       - wave.gif
%       - wave.avi



% sizes to match your screenshot: [Ny, Nx, Nt] = [501, 408, 1201]
Ny = 501; Nx = 408; Nt = 200;

load surfMesh.mat


X = xMesh; Y = yMesh; %Z = surfaceData1200;

dt = 0.01; t = (0:Nt-1)*dt;

A = 15.0;                         % peak height
a = 100.0;                         % width
omega = 2*pi*0.5;                % breathing frequency (Hz)

surfData = zeros(Ny, Nx, Nt, 'double');
G0 = exp(-(X.^2 + Y.^2)/a^2);    % static Gaussian

for k = 1:Nt
    surfData(:,:,k) = A * sin(omega*t(k)) .* G0;  % goes above and below zero
end

h5create('wave.h5','/u',[Ny Nx Nt],'ChunkSize',[Ny Nx 1],'Deflate',4);
for k = 1:Nt
    h5write('wave.h5','/u', surfData(:,:,k), [1 1 k], [Ny Nx 1]);
end

%k = 1; surf(x,y, surfData(:,:,k)); axis tight; shading interp; colormap parula

%% Save
for k = 1:size(surfData,3)
    img = mat2gray(surfData(:,:,k));    % scale to [0,1]
    [A,map] = gray2ind(img,256);
    if k == 1
        imwrite(A,map,'wave.gif','gif','LoopCount',Inf,'DelayTime',0.03);
    else
        imwrite(A,map,'wave.gif','gif','WriteMode','append','DelayTime',0.03);
    end
end

%% Save2
v = VideoWriter('wave.avi');   % or 'wave.mp4' if you have MPEG-4 support
v.FrameRate = 30;              % frames per second
open(v);

for k = 1:size(surfData, 3)
    surf(surfData(:,:,k));
    shading interp
    axis tight
    view(2)                    % 2D top view, or use view(3) for 3D
    caxis([-1 1])              % fix color range
    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v);

end