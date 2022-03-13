function [] = saveToFigureDir(fileNameID,figurePanelDir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
        saveas(gcf,fullfile(figurePanelDir,sprintf('%s.tif',fileNameID)))
   
             print(fullfile(figurePanelDir,fileNameID),'-depsc')

