function [timeStamp]=getTimeStamp()
	timeStamp=strrep(datestr(datetime),' ','_');
	timeStamp=strrep(timeStamp,':','-');
