function New_Posi = CorrectChaOrderCA1(FileBase,Shank,Posi)

% Shank; shank number
% Posi; recording site number in a shank.
%
EC16_1 = {'ec016.100_121','ec016.228_240','ec016.267_278','ec016.281_298','ec016.390_405','ec016.425_437','ec016.444_459'};
EC16_2 = {'ec016.479_487','ec016.532_540'};
EC16_3 = {'ec016.491_508'};
EC16_4 = {'ec016.577_590'};
GOR_1 = {'2006-6-7','6-08-2006','2006-6-12','2006-6-13'};
PIN_1 = {'2005-11-04'};
VVP_1 = {'2006-4-9','2006-4-10'};
VVP_2 = {'2006-4-18'};
G01_1 = {'g01_maze04_MS.001_007','g01_maze05_MS.001_002','g01_maze11_MS.001_004'};
KM01_1 = {'km01.004_011','km01.012_020'};
i01_1 = {'i01_maze01_MS.001_003','i01_maze02_MS.001_003','i01_maze04_MS.001_003','i01_maze05_MS.001_003','i01_maze15_MS.001_004'};
i01_2 = {'i01_maze15_MS.001_004'};

New_Posi = Posi;
if strcmp(FileBase(1:5),'ec013')
    if strcmp(FileBase,'ec013.451_470')
        if Shank==6
            New_Posi(Posi==4) = 5;
            New_Posi(Posi==5) = 4;
        elseif Shank==7
            New_Posi(Posi==6) = 7;
            New_Posi(Posi==7) = 8;
        end
    else
        if Shank==7
            New_Posi(Posi==4) = 5;
            New_Posi(Posi==5) = 4;
        elseif Shank==8
            New_Posi(Posi==6) = 7;
            New_Posi(Posi==7) = 8;
        end
    end
end

if strcmp(FileBase(1:5),'ec016')
    %if ismember(FileBase,EC16_1)
    if Shank==3
        New_Posi(Posi==7) = 8;
    elseif Shank==4
        New_Posi(Posi>1) =  New_Posi(Posi>1)+1;
    elseif Shank==9
        New_Posi(Posi==4) = 5;
        New_Posi(Posi==5) = 6;
        New_Posi(Posi==6) = 8;
    elseif Shank==10
        New_Posi(Posi==2) = 3;
        New_Posi(Posi==3) = 5;
        New_Posi(Posi==4) = 8;
    end

    %if ismember(FileBase,EC16_1)
        %if Shank==2
        %end
    if ismember(FileBase,EC16_2)
        if Shank==1
            New_Posi(Posi==7) = 8;
        end
    elseif ismember(FileBase,EC16_3)
        if Shank==1
            New_Posi(Posi==7) = 8;
        elseif Shank==7
            New_Posi(Posi>1) = New_Posi(Posi>1)+1;
        end
        %elseif Shank==10
        %New_Posi(Posi==2) = 3;
        %New_Posi(Posi==3) = 5;
        %end
        %else
        %error('aho')
    end
end

if ismember(FileBase,GOR_1)
    if Shank==9
        New_Posi(Posi==3) = 5;
        New_Posi(Posi==5) = 3;
    end
    New_Posi = 9-New_Posi;
end

if ismember(FileBase,G01_1)
    if Shank==4
        New_Posi(Posi==3) = 4;
        New_Posi(Posi==4) = 6;
        New_Posi(Posi==5) = 7;
        New_Posi(Posi==6) = 8;
        New_Posi(Posi==7) = 3;
        New_Posi(Posi==8) = 5;
    end
    if Shank==9
        New_Posi(Posi==4) = 8;
        New_Posi(Posi==8) = 4;
    end
    if Shank==16
        New_Posi(Posi==3) = 6;
    end
end

if ismember(FileBase,PIN_1)
    New_Posi = 9-New_Posi;
end

if ismember(FileBase,KM01_1)
    New_Posi = 9-New_Posi;
end

if ismember(FileBase,VVP_1)
    if Shank==2
        New_Posi(Posi==5) = 6;
        New_Posi(Posi==6) = 5;
    end
end

if ismember(FileBase,i01_1)
    if Shank==11 | Shank==13
        New_Posi = Posi.*2;
    end
    if Shank==7
        New_Posi(Posi>=4) = Posi(Posi>=4)+1;
    end
end
if ismember(FileBase,i01_2)
    if Shank==9
        New_Posi(Posi>=5) = Posi(Posi>=5)+1;
    end
end
        

%if ismember(FileBase,VVP_2)
    %if Shank==6
        %New_Posi(Posi==5) = 6;
        %New_Posi(Posi==6) = 5;
    %end
%end





