function [height, width, ratio] = calibrate(path, xrange, yrange)
    assert( exist(path, 'file') == 2, "video doesn't exist");
    
    v = VideoReader(path);
    frame = readFrame(v);
    frame = rgb2gray(frame);
    [height, width] = size(frame);
    
    % crop & binarize image
    frame_crop = frame(xrange(1):xrange(2), yrange(1):yrange(2));
    frame_bin = imbinarize(frame_crop, 0.5);
    frame_lbl = bwlabel(frame_bin, 4);
    lbls = unique(frame_lbl);
    lbls = lbls(2:length(lbls));
    n_lbls = length(lbls);
    
    % calculate mass center of each "spot"
    mc = zeros(n_lbls, 2, 'uint64');
    for i=1:length(lbls)
        [r, c] = find(frame_lbl == lbls(i));
        mc(i, :) = mean([r, c]);
    end
    mc = double(mc);
    
    % calculate # pixels between each adjacent "spot"
    % set unit-pixel-length as their avg.
    diffs = mc(2:n_lbls, :) - mc(1:n_lbls-1, :);
    unit_px = mean(hypot(diffs(:, 1), diffs(:, 2)));
    unit_cm = 1;
    
    % calculate cm - pixel conversion ratio
    ratio = unit_cm / unit_px;
    
end

