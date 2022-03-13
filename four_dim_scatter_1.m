% +------------------------------------------------------+
% |                4D Data Visualization                 |
% |              with MATLAB Implementation              | 
% |                                                      |
% | Author: M.Sc. Eng. Hristo Zhivomirov        02/16/13 | 
% +------------------------------------------------------+

clear, clc, close all

% form the axes
x = 0:0.5:10;                   % first dimension independent variable
y = 0:0.5:10;                   % second dimension independent variable
z = 0:0.5:10;                   % third dimension independent variable
[X, Y, Z] = meshgrid(x, y, z);  % form the 3D data grid

% form the data matrix - it is the fourth dimension
% the data could be imported from a file or could be generated via equation 
data = abs(cos(X) + cos(Y) + cos(Z));   

% plot the data
figure(1)
scatter3(X(:), Y(:), Z(:), 20, data(:), 'filled')
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
xlabel('X')
ylabel('Y')
zlabel('Z')
colorbar
