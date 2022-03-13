metaDataDir='./metadata';

metaDataFilePaths=getFilePathsRegex(metaDataDir,'*.mat');

for i=1:length(metaDataFilePaths)
    metaData=load(metaDataFilePaths{i});
    ratName=metaData.AnimalMetadata.AnimalName;
    ratTreatment=metaData.AnimalMetadata.Animal.ExperimentalTreatment;
    
    %ratTreatment=[ratTreatment '_' metaData.AnimalMetadata.Animal.GeneticLine];
  if(~strcmp(ratTreatment,'Sac') && ~strcmp(ratTreatment,'sac') )
    disp(sprintf('%s is %s',ratName,ratTreatment))
  end
end
