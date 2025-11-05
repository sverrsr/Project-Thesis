folder = 'screens';  % folder where .mat files are
files = dir(fullfile(folder, 'screen_*.mat'));
files = {files.name};
files = sort(files);

figure;
for k = 1:numel(files)
    data = load(fullfile(folder, files{k}));
    
    % Adjust if needed depending on what's inside:
    if isfield(data, 'screen_image')
        img = data.screen_image;
    elseif isfield(data, 'screen')
        img = data.screen.image;
    else
        fns = fieldnames(data);
        img = data.(fns{1});  % fallback if unsure
    end

    imagesc(img);
    axis image;
    colormap hot;
    colorbar;
    set(gca, 'YDir', 'normal');
    title(sprintf('Frame %d / %d', k, numel(files)));
    drawnow;
end
