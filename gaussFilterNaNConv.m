filtWidth = 7;
filtSigma = 5;
imageFilter=fspecial('gaussian',filtWidth,filtSigma);
dataFiltered = nanconv(data,imageFilter, 'nanout');