function [K,H,Pmax,Pmin] = curvature(X,Y,Z)
% Applies surfature function and make figures for exploratory data analysis


[K, H, Pmax, Pmin] = surfature(X, Y, Z);


% ---- Plot all results as subplots ----
figure;

% 1) Height map
subplot(2,2,1);
imagesc(X(1,:), Y(:,1), Z);
axis image;
set(gca, 'YDir', 'normal');
colormap(gca, turbo);
colorbar;
title('Surface elevation');
xlabel('X'); ylabel('Y');

% 2) Mean curvature (surfature)
subplot(2,2,2);
imagesc(X(1,:), Y(:,1), H);
axis image;
set(gca, 'YDir', 'normal');
colormap(gca, jet);
colorbar;
title('Mean curvature H');
xlabel('X'); ylabel('Y');

% 2) Mean curvature (surfature)
subplot(2,2,3);
imagesc(X(1,:), Y(:,1), K);
axis image;
set(gca, 'YDir', 'normal');
colormap(gca, jet);
colorbar;
title('Gaussian Curvature K');
xlabel('X'); ylabel('Y');

% 3) Histogram
subplot(2,2,4);
histogram(H(:), 100);
axis square;  % Make the axes square
title('Histogram of H');
xlabel('H value'); ylabel('Count');

% Find correlation of H and surface elevation
R_H = corr(H(:), Z(:), 'rows', 'complete');

end
