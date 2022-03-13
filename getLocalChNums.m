function [localChanNums]=getLocalChNums(cellPropFile,electrodeDiameter)
                localChNum=getChFromPath(cellPropFile);
                [row,col]=getSpatialRowColFromCh(cellPropFile,localChNum);

                electrodeRad=floor(electrodeDiameter/2);

                localChanNums=[];
                for localRow=(row-electrodeRad):(row+electrodeRad)
                        for localCol=(col-electrodeRad):(col+electrodeRad)
                                newCh=getChFromSpatialRowCol(cellPropFile,localRow,localCol);
                                localChanNums=[localChanNums newCh];
                        end
                end
