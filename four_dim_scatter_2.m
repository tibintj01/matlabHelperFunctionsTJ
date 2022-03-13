% +------------------------------------------------------+
% |                4D Data Visualization                 |
% |              with MATLAB Implementation              | 
% |                                                      |
% | Author: M.Sc. Eng. Hristo Zhivomirov        02/16/13 | 
% +------------------------------------------------------+

clear, clc, close all

% form the axes
x = 1:1:3;                      % first dimension independent variable
y = 1:1:3;                      % second dimension independent variable
z = 1:1:3;                      % third dimension independent variable
[X, Y, Z] = meshgrid(x, y, z);  % form the 3D data grid

% form the data matrix - it is the fourth dimension
% the data could be imported from a file or could be generated via equation 
data(:, :, 1) = [1  2  3;  4  5  6;  7  8  9];      % first page
data(:, :, 2) = [10 11 12; 13 14 15; 16 17 18];     % second page
data(:, :, 3) = [19 20 21; 22 23 24; 25 26 27];     % third page

% plot the data
figure(1)
scatter3(X(:), Y(:), Z(:), 500, data(:), 'filled')
grid3 on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
xlabel('X')
ylabel('Y')
zlabel('Z')
colorbar