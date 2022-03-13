function [sessName]=getSeizureStartTimes(cellPropFileName)


        sessName='error';

        if(length(findstr(cellPropFileName,'MG49_seizure36'))>0)
                %sessName='MG49_seizure36';
                sessName='MG49-seizure36';
        end
        if(length(findstr(cellPropFileName,'MG49-seizure43'))>0)
                sessName='MG49-seizure43';
        end
        if(length(findstr(cellPropFileName,'MG49-seizure45'))>0)
                sessName='MG49-seizure45';
        end
        if(length(findstr(cellPropFileName,'MG63_seizure1-4'))>0)
                %sessName='MG63_seizure1-4';
                sessName='MG63-seizure1-4';
        end
        if(length(findstr(cellPropFileName,'BW9_seizure1-3'))>0)
                %sessName='BW9_seizure1-3';
                sessName='BW9-seizure1-3';
        end
        if(length(findstr(cellPropFileName,'RIHE1'))>0)
                sessName='RIHE1';
        end
