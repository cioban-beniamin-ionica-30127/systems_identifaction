% Find optimal pair of lines to fit noisy data, one line on left side and one line on right side.  Separation/crossing x value is identified.
clc;	% Clear command window.
clear;	% Delete all variables.
close all;	% Close all figure windows except those created by imtool.
workspace;	% Make sure the workspace panel is showing.
fontSize = 18;
markerSize = 20;
%============================================================================================================
% FIRST CREATE SAMPLE X-Y DATA.  Skip (delete) this part if you already have your x and y data vectors.
% Define number of points.
numPoints = 100;
midIndex = round(numPoints/2);
% Load sample data.
s=load('product_16.mat')
x = s.time;
yNoisy = s.y;
% Done defining sample data.
%============================================================================================================
% Plot lines.
subplot(2, 2, 1);
plot(x, yNoisy, 'b.', 'MarkerSize', markerSize);
grid on;
xlabel('x', 'FontSize', fontSize);
ylabel('y', 'FontSize', fontSize);
title('Noisy Y vs. X data', 'FontSize', fontSize);
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
%================================================================================================================================================
% NOW SCAN ACROSS DATA FITTING LEFT AND RIGHT PORTION TO LINES.
% Assume the crossing point will be somewhere in the middle half of the points.
% If you go much more outward than that there are too few points to get a good line and the slopes of the lines will vary tremendously.
% Fit a line through the right and left parts and get the slopes.
% Keep the point where the slope difference is greatest.
numPoints = length(x);
index1 = round(0.25 * numPoints); % 25% of the way through.
index2 = round(0.75 * numPoints); % 75% of the way through.
% In other words, assume that we need at least 25 percent of the points to make a good estimate of the line.
% Obviously if we took only 2 or 3 points, then the slope could vary quite dramatically,
% so let's use at least 25% of the points to make sure we don't get crazy slopes.
% Initialize structure array
for k = 1 : numPoints
    lineData(k).slopeDifferences = 0;
    lineData(k).line1 = [0,0];
    lineData(k).line2 = [0,0];
end
for k = index1 : index2
    % Get data in left side.
    x1 = x(1:k);
    y1 = yNoisy(1:k);
    % Fit a line through the left side.
    coefficients1 = polyfit(x1, y1, 1); % The slope is coefficients1(1).
    % Get data in right side.
    x2 = x(k+1:end);
    y2 = yNoisy(k+1:end);
    % Fit a line through the left side.
    coefficients2 = polyfit(x2, y2, 1); % The slope is coefficients2(1).
    
    % Compute difference in slopes, and store in structure array along with line equation coefficients.
    lineData(k).slopeDifferences = abs(coefficients1(1) - coefficients2(1));
    lineData(k).line1 = coefficients1;
    lineData(k).line2 = coefficients2;
end
% Find index for which slope difference is greatest.
slopeDifferences = [lineData.slopeDifferences]; % Extract from structure array into double vector of slope differences only
% slope1s = struct2table(lineData.line1); % Extract from structure array into double vector of slopes only
% slope2s = [lineData.line2(1)]; % Extract from structure array into double vector of slopes only
[maxSlopeDiff, indexOfMaxSlopeDiff] = max(slopeDifferences)
% Plot slope differences.
subplot(2, 2, 2);
plot(slopeDifferences, 'b.', 'MarkerSize', markerSize);
xlabel('Index', 'FontSize', fontSize);
ylabel('Slope', 'FontSize', fontSize);
grid on;
caption = sprintf('Slope Differences Maximum at Index = %d, x value = %.2f', indexOfMaxSlopeDiff, x(indexOfMaxSlopeDiff));
title(caption, 'FontSize', fontSize);
% Mark it with a red line.
line([indexOfMaxSlopeDiff, indexOfMaxSlopeDiff], [0, maxSlopeDiff], 'Color', 'r', 'LineWidth', 2);
% Show everything together all on one plot.
% Plot lines.
subplot(2, 2, 3:4);
plot(x, yNoisy, 'b.', 'MarkerSize', markerSize);
grid on;
xlabel('x', 'FontSize', fontSize);
ylabel('y', 'FontSize', fontSize);
hold on;
% Use the equation of line1 to get fitted/regressed y1 values.
slope1 = lineData(indexOfMaxSlopeDiff).line1(1);
intercept1 = lineData(indexOfMaxSlopeDiff).line1(2);
y1Fitted = slope1 * x + intercept1;
% Plot line 1 over/through data.
plot(x, y1Fitted, 'r-', 'LineWidth', 2);
% Use the equation of line2 to get fitted/regressed y2 values.
slope2 = lineData(indexOfMaxSlopeDiff).line2(1);
intercept2 = lineData(indexOfMaxSlopeDiff).line2(2);
y2Fitted = slope2 * x + intercept2;
% Plot line 2 over/through data.
plot(x, y2Fitted, 'r-', 'LineWidth', 2);
%================================================================================================================================================
% FIND THE CROSSING POINT.  IT IS WHERE THE Y VALUES ARE EQUAL.
% Just set the equations equal to each other and solve for x.  
% So y1Fitted = y2Fitted, which means (slope1 * xc + intercept1) = (slope2 * xc + intercept2).  Solving for xc gives:
xc = (intercept2 - intercept1) / (slope1 - slope2);
y1c = slope1 * xc + intercept1; % Will be the same as y2c.
y2c = slope2 * xc + intercept2; % Will be the same as y1c.
% Mark crossing with a magenta line.
% Draw a line up from the x axis to the crossing point.
line([xc, xc], [0, y1c], 'Color', 'm', 'LineWidth', 2);
% Draw a line over from the y axis to the crossing point.
line([0, xc], [y1c, y1c], 'Color', 'm', 'LineWidth', 1);
caption = sprintf('Data with left and right lines overlaid.  Lines cross at (x,y) = (%.4f, %.4f)', xc, y1c);
title(caption, 'FontSize', fontSize);
message1 = sprintf('Left  Equation: y = %.3f * x + %.3f', slope1, intercept1);
message2 = sprintf('Right Equation: y = %.3f * x + %.3f', slope2, intercept2);
message = sprintf('%s\n%s', message1, message2);
fprintf('%s\n', message);
text(5, 100, message, 'Color', 'r', 'FontSize', 15, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
uiwait(helpdlg(message));