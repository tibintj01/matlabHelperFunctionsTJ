function newString = removeUnderscores(oldString)

underscoreInd = strfind(oldString, '_');
oldString(underscoreInd) = ' ';
newString = oldString;