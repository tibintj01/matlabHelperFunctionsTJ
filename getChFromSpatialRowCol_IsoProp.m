function [chNum]=getChFromSpatialRowCol_IsoProp(isoProps,row,col)
	 	
		session_list_human_seizures_NEW_FOR_ELLEN_SORTED_DATA_2017_TJ

		szSubject=dataInfo(isoProps.sessionIdx).subject;

		patientID=szSubjectToPatientID(szSubject);

	if(strcmp(patientID,'BW9'))
                patientID='MG49'; %no BW9.xls but all same channel mapping (Omar communication, Jan 20, 2019)
        end
                
		[arrayMap,electrodeXY,electrodeImp]=neuroportArrayData(patientID);

		chNum=-1;
		for ch=1:96
			if(electrodeXY(ch,1)==col && ((10-electrodeXY(ch,2))+1 ==row))
				chNum=ch;
			end
		end
