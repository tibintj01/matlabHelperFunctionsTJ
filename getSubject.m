function [sessName]=getSeizureStartTimes(cellPropFileName)


        sessName='error';

        if(length(findstr(cellPropFileName,'MG49'))>0)
                %sessName='MG49_seizure36';
                sessName='MG49';
        end
        if(length(findstr(cellPropFileName,'MG63'))>0)
                %sessName='MG63_seizure1-4';
                sessName='MG63';
        end
        if(length(findstr(cellPropFileName,'BW9'))>0)
                %sessName='BW9_seizure1-3';
                sessName='BW9';
        end
        if(length(findstr(cellPropFileName,'RIHE1'))>0)
                sessName='RIHE1';
        end
