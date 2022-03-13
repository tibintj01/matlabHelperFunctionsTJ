function [row,col]=getSpatialRowColFromIsoProps(isoProps,ch)
	session_list_human_seizures_NEW_FOR_ELLEN_SORTED_DATA_2017_TJ

        szSubject=dataInfo(isoProps.sessionIdx).subject;

	patientID=szSubjectToPatientID(szSubject);

	if(strcmp(patientID,'BW9'))
		patientID='MG49'; %no BW9.xls but all same channel mapping (Omar communication, Jan 20, 2019)
	end	
	[arrayMap,electrodeXY,electrodeImp]=neuroportArrayData(patientID);

	%electrodeXY
	%row=electrodeXY(ch,1);
	%col=electrodeXY(ch,2);
	
	col=electrodeXY(ch,1);
	row=(10-electrodeXY(ch,2))+1;
end
