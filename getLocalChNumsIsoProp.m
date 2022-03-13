function [localChanNums]=getLocalChNumsIsoProp(isoProps,electrodeDiameter)
		session_list_human_seizures_NEW_FOR_ELLEN_SORTED_DATA_2017_TJ

                localChNum=isoProps.ch;

                [row,col]=getSpatialRowColFromIsoProps(isoProps,localChNum);

                electrodeRad=floor(electrodeDiameter/2);

                localChanNums=[];
                for localRow=(row-electrodeRad):(row+electrodeRad)
                        for localCol=(col-electrodeRad):(col+electrodeRad)
                                newCh=getChFromSpatialRowCol_IsoProp(isoProps,localRow,localCol);
                                localChanNums=[localChanNums newCh];
                        end
                end
