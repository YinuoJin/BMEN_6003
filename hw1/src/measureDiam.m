% Load frames from given video, select the appropriate frame to calculate 
% EDV & ESV of the ventricular cavity diameter in the given view
% Modified from the following tutorial / scripts:
% i). Interactive distance measure: 
% https://www.mathworks.com/help/images/measure-distances-in-images.html
% ii). Display sequence of images:
% https://www.mathworks.com/matlabcentral/answers/195418-how-to-show-a-sequence-of-images


function [diam_ed, diam_es] =  measureDiam(path, height, width, scale)
    assert( exist(path, 'file') == 2, "video doesn't exist");

    % (1). Variable initialization
    global pos_ed;  % end-to-end mouse pixels for EDV measure positions
    global pos_es;  % end-to-end mouse pixels for ESV measure positions
    
    pos_ed = [];
    pos_es =[];
    
    % (2). Load sample frame & get specs  
    v = VideoReader(path);
    I = zeros(height, width, v.NumFrames, 'uint8');

    my_data.Units = 'cm';
    my_data.MaxValue = hypot(height, width);
    my_data.Colormap = hsv;
    my_data.ScaleFactor = scale;

    % (3). Load frames & measurement
    i = 1;
    while hasFrame(v)
        frame = readFrame(v);
        frame = rgb2gray(frame);
        I(:, :, i) = frame;
        i = i+1;
    end

    count = 1;
    
    h_img = image(I(:, :, count));
    h_img.ButtonDownFcn = @(~,~) draw(h_img.Parent,my_data);

    hb1 = uicontrol('Style', 'PushButton', 'String', 'Next', ...
      'Callback', @nextFrame);

    hb2 = uicontrol('Style', 'PushButton', 'String', 'Prev', ...
      'Callback', @prevFrame);

    hb1.Position(2) = hb2.Position(2)+1.1*hb2.Position(4);

    while isempty(pos_ed) || isempty(pos_es)
        pause(.1);
    end

    diam_ed = calcDiam(pos_ed, my_data.ScaleFactor);
    diam_es = calcDiam(pos_es, my_data.ScaleFactor);
    
    
    function nextFrame(src, evt)
        if count < size(I, 3)
            count = count + 1;
            h_img.CData = I(:, :, count);
        end
    end


    function prevFrame(src, evt)
        if count > 1
            count = count - 1;
            h_img.CData = I(:, :, count);
        end
    end

end


% Helper functions
%%%%%%%%%%%%%%%%%%%%%%%

function draw(hAx,my_data)
    global pos_ed;
    global pos_es;

    % Create a line ROI object. Specify the initial color of the line and
    % store the |my_data| structure in the |UserData| property of the ROI.
    h = images.roi.Line('Color',[0, 1, 0],'UserData',my_data);

    % Set up a listener for movement of the line ROI.
    addlistener(h,'MovingROI',@updateLabel);

    % Set up a listener for clicks on the line ROI.
    addlistener(h,'ROIClicked',@updateUnits);

    % Get the current mouse location from the |CurrentPoint| property of the
    % axes and extract the _x_ and _y_ coordinates.
    cp = hAx.CurrentPoint;
    cp = [cp(1,1) cp(1,2)];

    % Begin drawing the ROI from the current mouse location. Using the
    % |beginDrawingFromPoint| method, you can draw multiple ROIs.

    hROIs = findobj(hAx,'Type','images.roi.Line');
    if numel(hROIs) < 2  
        h.beginDrawingFromPoint(cp);
    else
        % LIFO order
        pos_ed = get(hROIs(2), 'Position');
        pos_es = get(hROIs(1), 'Position');
        set(hROIs, 'UserData', pos_ed);
        set(hROIs, 'UserData', pos_es);
    end

end


function updateLabel(src, evt)
    % Get the current line position.
    pos = evt.Source.Position;

    % Determine the length of the line.
    diffPos = diff(pos);
    mag = hypot(diffPos(1),diffPos(2));

    % Choose a color from the colormap based on the length of the line. The
    % line changes color as it gets longer or shorter.
    color = src.UserData.Colormap(ceil(256*(mag/src.UserData.MaxValue)),:);

    % Apply the scale factor to line length to calibrate the measurements.
    mag = mag*src.UserData.ScaleFactor;

    % Update the label.
    set(src,'Label',[num2str(mag,'%30.1f') ' ' src.UserData.Units],'Color',color);

end


function updateUnits(src, evt)
    % When you double-click the ROI label, the example opens a popup dialog box
    % to get information about the actual distance. Use this information to
    % scale all line ROI measurements.
    if strcmp(evt.SelectionType,'double') && strcmp(evt.SelectedPart,'label')

        % Display the popup dialog box.
        answer = inputdlg({'Known distance','Distance units'},...
            'Specify known distance',[1 20],{'10','meters'});

        % Determine the scale factor based on the inputs.
        num = str2double(answer{1});

        % Get the length of the current line ROI.
        pos = src.Position;
        diffPos = diff(pos);
        mag = hypot(diffPos(1),diffPos(2));

        % Calculate the scale factor by dividing the known length value by the
        % current length, measured in pixels.
        scale = num/mag;

        % Store the scale factor and the units information in the |my_data|
        % structure.
        my_data.Units = answer{2};
        my_data.MaxValue = src.UserData.MaxValue;
        my_data.Colormap = src.UserData.Colormap;
        my_data.ScaleFactor = scale;

        % Reset the data stored in the |UserData| property of all existing line
        % ROI objects. Use |findobj| to find all line ROI objects in the axes.
        hAx = src.Parent;
        hROIs = findobj(hAx,'Type','images.roi.Line');
        set(hROIs,'UserData',my_data);

        % Update the label in each line ROI object, based on the information
        % collected in the input dialog.
        for i = 1:numel(hROIs)
            pos = hROIs(i).Position;
            diffPos = diff(pos);
            mag = hypot(diffPos(1),diffPos(2));
            set(hROIs(i),'Label',[num2str(mag*scale,'%30.1f') ' ' answer{2}]);

        end

        % Reset the |ButtonDownFcn| callback function with the current |my_data|
        % value.
        h_img = findobj(hAx,'Type','image');
        h_img.ButtonDownFcn = @(~,~) draw(hAx,my_data);

    end

end



function diam = calcDiam(pos, scale)
    diffs = diff(pos);
    diam = scale*hypot(diffs(1), diffs(2));
end
