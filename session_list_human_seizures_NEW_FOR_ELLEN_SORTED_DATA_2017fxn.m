function [dataInfo]=session_list_human_seizures_NEW_FOR_ELLEN_SORTED_DATA_2017fxn()

%% SESSION LIST TO ANALYZE HUMAN SEIZURES
%  Use this for all paper related code
%  The key differences from "session_list_human_seizures" is that this file works on one ns5 at at time, even when rethresholding

%% BLACKROCK DATA TYPE
%  The choice may influence exactly how the data is processed (e.g. how to
%  construct array maps of the data).
NPORT = 1; MWIRE = 2; MGRID = 3; LAMIN = 4;

tjStartAdj_BW9_Sz1=-15;
tjEndAdj_BW9_Sz1=-1;

tjStartAdj_BW9_Sz2=12;
tjEndAdj_BW9_Sz2=-11;

tjStartAdj_BW9_Sz3=5;
tjEndAdj_BW9_Sz3=3;

tjStartAdj_MG49_43=1;
tjEndAdj_MG49_43=16;

tjEndAdj_RIHE1_Sz1=133;

changeMe = 0;
clear dataInfo;

%% SEIZURE WAVE BY WAVE ANALYSIS

%% SEIZURE SESSION LIST
%  1.  MG49 Sz 36, USE THIS (ELLEN SORTED)
%  2.  MG49 Sz 43, USE THIS (ELLEN SORTED)
%  3.  MG49 Sz 45, USE THIS (ELLEN SORTED)
%  4.  MG63  Sz 1-4, concat, rethresh
%  5.  BW9 Sz 1-3, USE THIS (ELLEN SORTED)
%  6.  RIHE1 Sz 2 + 3 USE THIS (JASON SORTED)
%  7.  MG29 Sz 22, DO NOT USE THIS - DIFFERENT KIND OF SEIZURE

%changed 43b quality from 1 to 5- Tibin, 9/1/18
%% id = 1 MG49 SZ36 - ELLEN HAS SORTED NS5
changeMe = changeMe + 1;
    %'subject', 'MG49', ...
dataInfo(changeMe) = struct( ...
    'subject', 'MG49_seizure36', ...
    'cellsFilename', '', ... % Include suffix. Sorted Units. the name of the NEX (NEUROEXPLORER) file that contains the final clusters to use for this dataset
    'descriptiveFilename', 'MG49_seizure36_20110612-151513-010_011_012', ... % This is just a shortform for the the cellsFilename when saving heavily processed data
    'rethresholded', 1, ... % 0=clustered as original NEX, 1=rethresholded (with or without concatenation... this has 1 channel per file
    'cellChannel',  [3	4	4	5	6	7	7	8	8	9	12	15	16	16	17	18	18	19	20	20	21	21	22	23	23	24	25	26	27	28	29	30	31	31	32	33	34	35	36	38	39	40	41	42	43	43	43	44	44	45	46	47	48	49	50	51	51	52	53	53	55	56	58	59	60	64	70	84	84	86	88	91	91	92	92	93	94	95	96], ...    
    'cellIsInterneuron', [], ...  [ 0  1  0], ...
    'cellQuality',       [1	4	4	2	1	2	3	4	3	2	3	3	4	5	3	4	3	2	4	3	3	3	4	3	1	4	4	4	4	4	3	5	2	2	3	2	1	1	2	3	3	4	3	3	5	5	4	2	3	1	3	4	3	1	1	4	4	3	3	1	3	3	1	1	4	2	5	5	4	3	4	3	1	5	4	4	4	3	4], ...  [ 4  4  3], ...
    ... % Seizure Info:
    'numSz', 1, ...
    'szStartTime', 3120, ...
    'szEndTime', 3200, ...
    'szPlotStartTime', 3100, ...
    'szPlotEndTime', 3220, ...
    ... % Matching NS5 Info:
    'numNS5Files', 3, ... % the number of LFP files that are to be analyzed as one
    'ns5Filenames', {{'20110612-151513-010', '20110612-151513-011', '20110612-151513-012'}}, ... % a cell array
    'ns5VisChannels', [1:96], ... % As seen on Blackrock Central. Leave blank to use all
    'ns5SeqChannels', [], ... % Sequential channels, 
    'ns5TrigVisChannel', 129, ... % The trigger channel, as seen on Blackrock Central
    'ns5TrigSeqChannel', 97, ...  % The sequenctial number for the trigger channel
    ... % Matching ECoG Info:
    'numEDFFiles', 3, ... % set to 0 if no corresponding ECoG files, or bad triggers, or BrainGate data. EcoG filenames: one entry for each NS5 that it matches
    'edfFilenames', {{'MG49_Seizure36', 'MG49_Seizure36', 'MG49_Seizure36'}}, ... 
    'edfClosestVisChannels', {{'GR35'}}, ...
    'edfChannels', [1:21 23:32 34:83 85:112], ... % GR32 is noisy. GR21 (seq22), PSTS3 (seq84) are missing. 2-65 = grid, sequential channels to process
    'edfTrigVisChannel', 'Event', ...
    'edfTrigSeqChannel', 1 ...
);

% Below is an entry for the NEX version of seizure 36. We sorted this for
% comparison sake and including this note in case you want to look at this
% data for any reason to compare concat + rethresholded sorting vs NEX
% sorting
% %% id = 2 MG49 SZ36 Cell 51a - ELLEN HAS SORTED THIS NEX FILE, but will not be used (will use above NS5 instead)
% changeMe = changeMe + 1;
% dataInfo(changeMe) = struct( ...
%     'subject', 'MG49', ...
%     'cellsFilename', '20110612-151513-011.nex', ... % Include suffix. Sorted Units. the name of the NEX (NEUROEXPLORER) file that contains the final clusters to use for this dataset
%     'descriptiveFilename', 'MG49_seizure36', ... % This is just a shortform for the the cellsFilename when saving heavily processed data
%     'rethresholded', 1, ... % 0=clustered as original NEX, 1=rethresholded (with or without concatenation... this has 1 channel per file
%     'cellChannel',       [], ... [43 43 43], ...
%     'cellIsInterneuron', [], ...  [ 0  1  0], ...
%     'cellQuality',       [], ...  [ 4  4  3], ...
%     ... % Seizure Info:
%     'numSz', 1, ...
%     'szStartTime', 720, ...
%     'szEndTime', 800, ...
%     'szPlotStartTime', 700, ...
%     'szPlotEndTime', 820, ...
%     ... % Matching NS5 Info:
%     'numNS5Files', 1, ... % the number of LFP files that are to be analyzed as one
%     'ns5Filenames', {{'20110612-151513-011'}}, ... % a cell array
%     'ns5VisChannels', [1:96], ...[8 13 16 30 40 43 51], ... % As seen on Blackrock Central. Leave blank to use all
%     'ns5SeqChannels', [], ... % Sequential channels, 
%     'ns5TrigVisChannel', 129, ... % The trigger channel, as seen on Blackrock Central
%     'ns5TrigSeqChannel', 97, ...  % The sequenctial number for the trigger channel
%     ... % Matching ECoG Info:
%     'numEDFFiles', 3, ... % set to 0 if no corresponding ECoG files, or bad triggers, or BrainGate data. EcoG filenames: one entry for each NS5 that it matches
%     'edfFilenames', {{'', '', ''}}, ... 
%     'edfClosestVisChannels', {{'GR35'}}, ...
%     'edfChannels', [1:21 23:32 34:83 85:112], ... % GR32 is noisy. GR21 (seq22), PSTS3 (seq84) are missing. 2-65 = grid, sequential channels to process
%     'edfTrigVisChannel', 'Event', ...
%     'edfTrigSeqChannel', 1 ...
% );

%% id = 2 MG49 SZ43 all channels - CONCATENATED AND GIVEN TO ELLEN, SHE HAS FINISHED SORTING
    %'subject', 'MG49-seizure43', ...
changeMe = changeMe + 1;
dataInfo(changeMe) = struct( ...
    'subject', 'MG49_seizure43', ...
    'cellsFilename', '', ... % Include suffix. Sorted Units. the name of the NEX (NEUROEXPLORER) file that contains the final clusters to use for this dataset
    'descriptiveFilename', 'MG49_seizure43_20110612-151513-024_025_026', ... % This is just a shortform for the the cellsFilename when saving heavily processed data
    'rethresholded', 1, ... % 0=clustered s originl NEX, 1=rethresholded (with or without conctention... this hs 1 chnnel per file
    'cellChannel', [1	3	4	5	6	7	8	8	9	10	11	12	13	13	14	16	16	16	17	18	19	20	21	22	23	23	24	25	26	26	27	27	28	28	29	30	30	31	32	32	33	33	34	35	36	37	37	38	39	40	40	41	42	43	43	43	43	45	46	46	47	48	49	50	51	51	52	52	53	54	55	56	58	58	59	59	60	60	61	62	63	63	64	64	66	68	73	74	75	76	79	80	82	83	84	85	85	86	87	88	91	92	92	92	93	93	94	94	95	96], ...
    'cellIsInterneuron', [], ...
    'cellQuality', [2	3	3	3	1	3	3	3	3	2	3	2	4	4	4	5	5	3	3.5	3	2	3	3.5	3.5	3	2	4	3	4	5	3.5	1.5	4	5	3	5	5	3	4	4	1	2	3	1.5	3	1	2	3	3	3	3.5	5	3.5	5	4	2	3	2	3	1	5	3.5	4	4	3.5	3.5	3.5	3.5	3.5	4	3	4	3.5	4	4	3.5	4	4	3.5	4	3.5	2	4	3	1	1	2	2	2	1	2	2	5	3	5	4	3	3	2	3.5	3.5	3	4	2	5	4	4	3.5	2	3.5], ...
    ... % Seizure Info:
    'numSz', 1, ...
    'szStartTime', (560+2400+tjStartAdj_MG49_43), ...
    'szEndTime', (647+2400+tjEndAdj_MG49_43), ...
    'szPlotStartTime', (546+2400), ...
    'szPlotEndTime', (676+2400), ...
    ... % Matching NS5 Info:
    'numNS5Files', 3, ... % the number of LFP files that are to be analyzed as one
    'ns5Filenames', {{'20110612-151513-024', '20110612-151513-025', '20110612-151513-026'}}, ... % a cell array
    'ns5VisChannels', [1:96], ... % As seen on Blackrock Central. Leave blank to use all
    'ns5SeqChannels', [], ... % Sequential channels, 
    'ns5TrigVisChannel', 129, ... % The trigger channel, as seen on Blackrock Central
    'ns5TrigSeqChannel', 97, ...  % The sequenctial number for the trigger channel
    ... % Matching ECoG Info:
    'numEDFFiles', 3, ... % set to 0 if no corresponding ECoG files, or bad triggers, or BrainGate data. EcoG filenames: one entry for each NS5 that it matches
    'edfFilenames', {{'', '', ''}}, ... 
    'edfClosestVisChannels', {{'GR35'}}, ...
    'edfChannels', [1:21 23:32 34:83 85:112], ... % GR32 is noisy. GR21 (seq22), PSTS3 (seq84) are missing. 2-65 = grid, sequential channels to process
    'edfTrigVisChannel', 'Event', ...
    'edfTrigSeqChannel', 1 ...
);

%% id = 3 MG49 SZ45 all channels - CONCATENATED AND GIVEN TO ELLEN, SHE HAS FINISHED SORTING
changeMe = changeMe + 1;
dataInfo(changeMe) = struct( ...
    'subject', 'MG49-seizure45', ...
    'cellsFilename', '', ... % Include suffix. Sorted Units. the name of the NEX (NEUROEXPLORER) file that contains the final clusters to use for this dataset
    'descriptiveFilename', 'MG49_seizure45_20110613-112805-068_069_070', ... % This is just a shortform for the the cellsFilename when saving heavily processed data
    'rethresholded', 1, ... % 0=clustered as original NEX, 1=rethresholded (with or without concatenation... this has 1 channel per file
    'cellChannel',       [3     4       5       6       7       7       8       9       9       9       12      12      12      13      14      14      15      16      16      17      18      18      19      19      20      20      21      21      22      22      22      22      23      24      24      25      25      26      26      27      27      28      28      29      29      30      30      31      31      31      32      33      34      36      37      38      39      39      40      41      42      43      43      46      47      47      47      48      49      50      51      52      53      54      54      55      56      56      57      57      58      58      59      59      60      60      61      61      62      62      63      63      64      64      66      68      70      73      74      75      76      78      79      81      82      83      84      85      86      86      87      88      91      91      92      92      93      93      94      94      95      96      96      96], ...
    'cellIsInterneuron', [], ...
    'cellQuality',       [2     3       3       3       4       3       2       4       3       3       4       3       1.5     4       4.5     4       4       2       3.5     3       4       3       2       3       4       3       3.5     3       1       2       3       4       3       5       4       3.5     3       4       3       3       2       4       3.5     3       3.5     4       4       4       4       2       3       2       4.5     3       1       3       4       3       3       5       3.5     4       4       4       5       4       1.5     2.5     2       3.5     4       5       4       4       4       3       2       3       4       3       3.5     3.5     4       3       4       3       4       2.5     5       1       5       3.5     5       4       1       4       2       2       1       3       1.5     1       3.5     2       4       2       2       3       3       3       2       3       4       3       4.5     3.5     4       4       4       3.5     2       4       3       3], ...
    ... % Seizure Info:
    'numSz', 1, ...
    'szStartTime', (255+2400), ...
    'szEndTime', (355+2400), ...
    'szPlotStartTime', (235+2400), ...
    'szPlotEndTime', (375+2400), ...
    ... % Matching NS5 Info:
    'numNS5Files', 3, ... % the number of LFP files that are to be analyzed as one
    'ns5Filenames', {{'20110613-112805-068', '20110613-112805-069', '20110613-112805-070'}}, ... % a cell array
    'ns5VisChannels', [1:96], ... % As seen on Blackrock Central. Leave blank to use all
    'ns5SeqChannels', [], ... % Sequential channels, 
    'ns5TrigVisChannel', 129, ... % The trigger channel, as seen on Blackrock Central
    'ns5TrigSeqChannel', 97, ...  % The sequenctial number for the trigger channel
    ... % Matching ECoG Info:
    'numEDFFiles', 3, ... % set to 0 if no corresponding ECoG files, or bad triggers, or BrainGate data. EcoG filenames: one entry for each NS5 that it matches
    'edfFilenames', {{'', '', ''}}, ... 
    'edfClosestVisChannels', {{'GR35'}}, ...
    'edfChannels', [1:21 23:32 34:83 85:112], ... % GR32 is noisy. GR21 (seq22), PSTS3 (seq84) are missing. 2-65 = grid, sequential channels to process
    'edfTrigVisChannel', 'Event', ...
    'edfTrigSeqChannel', 1 ...
);


%% id = 4 MG63 SZ1-4 - KEEP AS A SEQUENCE OF 6 - ELLEN SORTED
%  The seizure times below fall in file 24 and file 27, so it makes sense to
%  have all 6 sequential files from 23 to 28
changeMe = changeMe + 1;
dataInfo(changeMe) = struct( ...
    'subject', 'MG63_seizure1-4', ...
    'cellsFilename', '', ... % Include suffix. Sorted Units. the name of the NEX (NEUROEXPLORER) file that contains the final clusters to use for this dataset
    'descriptiveFilename', 'MG63_seizures1-4', ... % This is just a shortform for the the cellsFilename when saving heavily processed data
    'rethresholded', 1, ... % 0=clustered s originl NEX, 1=rethresholded (with or without conctention... this hs 1 chnnel per file
    'cellChannel',       [1	1	2	2	3	4	4	4	5	5	6	7	8	8	8	10	10	10	11	12	12	13	14	15	17	18	18	19	19	20	20	21	23	23	25	27	28	30	31	34	35	35	36	36	37	37	37	38	38	39	39	39	40	40	41	42	42	43	43	44	44	45	45	45	46	46	46	47	47	48	48	49	49	50	50	51	52	52	53	53	53	53	54	54	55	55	56	56	57	57	58	58	58	59	60	61	62	62	63	63	65	66	68	68	69	69	70	70	71	71	72	72	74	74	75	75	76	76	77	77	78	78	79	79	80	81	81	82	82	83	83	84	84	85	85	87	87	88	88	88	90	90	91	91	92	94	94	95	95], ...
    'cellIsInterneuron', [], ...
    'cellQuality',       [5	3	3	4	1.5	4	3	2	3.5	2	3	3	4	3	2	5	4	2	2	1	1	1.5	2	2	3.5	2	3	2	3	4	3	3	1	2	2.5	1	1	2	1.5	2	3	3	3	2	3.5	2	3	4	2.5	3.5	3	1.5	3	3	3	4	3	3	4	4	2	4	3	1.5	5	2.5	2.5	4	2.5	4	2	2	4	4.5	2.5	4	4	3	3	3	2	2	3	4.5	5	3	4	2.5	5	3	5	3	1.5	3	4	2	4	3	4	2	3	4	4	2	4	3	4	2	3	3	4	3	4	4	3	4	4	4	4	3	3	1	3	1.5	2.5	1	3.5	3	1	4	1	4	2	4	2	4	2	4	4	2	4	2	4	2	1	5	1	4	2], ...  
    ... % Seizure Info:
    'numSz', 2, ... 
    'szStartTime', [6532 12984], ... % only 2 of the seizures were said to be full blown clinical seizures. These are their start times.
    'szEndTime', [6632 13084], ...
    'szPlotStartTime', [6512 12964], ...
    'szPlotEndTime', [6652 13104], ...
    ... % Matching NS5 Info:
    'numNS5Files', 6, ... % the number of LFP files that are to be analyzed as one
    'ns5Filenames', {{'05-103012-023', '05-103012-024', '05-103012-025', '05-103012-026', ...
                      '05-103012-027', '05-103012-028'}}, ... % a cell array
    'ns5VisChannels', 1:96, ... [17 38 91], ... % As seen on Blackrock Central. Leave blank to use all
    'ns5SeqChannels', [], ... % Sequential channels, 
    'ns5TrigVisChannel', 129, ... % The trigger channel, as seen on Blackrock Central
    'ns5TrigSeqChannel', 97, ...  % The sequenctial number for the trigger channel
    ... % Matching ECoG Info:
    'numEDFFiles', 3, ... % set to 0 if no corresponding ECoG files, or bad triggers, or BrainGate data. EcoG filenames: one entry for each NS5 that it matches
    'edfFilenames', {{'', '', ''}}, ... 
    'edfClosestVisChannels', {{''}}, ...
    'edfChannels', [], ... 
    'edfTrigVisChannel', 'Event', ...
    'edfTrigSeqChannel', 1 ...
);

%% id = 5 BW9 SZ1-3 CONCATENATED, ELLEN SORTED, cellChannel and cellQuality filled in below by Omar based on Ellen's excel file notes
changeMe = changeMe + 1;
dataInfo(changeMe) = struct( ...
    'subject', 'BW9_seizure1-3', ...
    'cellsFilename', '', ... % Include suffix. Sorted Units. the name of the NEX (NEUROEXPLORER) file that contains the final clusters to use for this dataset
    'descriptiveFilename', 'BW9_seizures1-3', ... % This is just a shortform for the the cellsFilename when saving heavily processed data
    'rethresholded', 1, ... % 0=clustered as original NEX, 1=rethresholded (with or without concatenation... this has 1 channel per file
    'cellChannel', [2	2	2	8	10	11	11	12	14	14	16	16	22	23	24	24	27	28	29	29	30	31	37	37	51	53	55	56	57	59	63	64	65	66	69	69	72	87	88	90	91	93	95], ...     
    'cellIsInterneuron', [], ...
    'cellQuality', [1	2	2	3	3	2	3	3	3	3	3	2	3	3	1	2	3	1.5	2	3	1	2.5	0.5	0.5	1	2.5	3	2	1	3	2	2	1	1	3	1	1	1.5	2	3	2	3.5	1],...
    ... % Seizure Info:
    'numSz', 3, ... % the number of seizures in this session
    'szStartTime', [1500+tjStartAdj_BW9_Sz1 3850+tjStartAdj_BW9_Sz2  6725+tjStartAdj_BW9_Sz3], ... % need to determine exact times still based on visualization of LFP data
    'szEndTime', [1580+tjEndAdj_BW9_Sz1 3925+tjEndAdj_BW9_Sz2 6785+tjEndAdj_BW9_Sz3], ... % need to determine exact times still based on visualization of LFP data
    'szPlotStartTime', 0, ...
    'szPlotEndTime', 0, ...
    ... % Matching NS5 Info:
    'numNS5Files', 3, ... % the number of LFP files that are to be analyzed as one
    'ns5Filenames', {{'20081229-180214-101', '20081229-180214-102', '20081229-180214-103'}}, ... % a cell array
    'ns5VisChannels', [1:96], ... % As seen on Blackrock Central. Leave blank to use all
    'ns5SeqChannels', [], ... % Sequential channels, 
    'ns5TrigVisChannel', 129, ... % The trigger channel, as seen on Blackrock Central
    'ns5TrigSeqChannel', 97, ...  % The sequenctial number for the trigger channel
    ... % Matching ECoG Info:
    'numEDFFiles', 3, ... % set to 0 if no corresponding ECoG files, or bad triggers, or BrainGate data. EcoG filenames: one entry for each NS5 that it matches
    'edfFilenames', {{'', '', ''}}, ... 
    'edfClosestVisChannels', {{''}}, ...
    'edfChannels', [], ...
    'edfTrigVisChannel', 'Event', ...
    'edfTrigSeqChannel', 1 ...
);

%JANUARY 23 2018 - SzEndTime said 852, switched to 732+133 based on heatmap/raster - Tibin
%% id = 6 RIHE1 Sz2 to Sz3 with awake or asleep (unconfirmed) in between. 11 Files
%  XXX NOISY INTERVAL FROM 9515-10081
%  SORTING STATUS: Jason Finished Sorting, USING ORIGINAL JASON SORTED -
%  NO ELLEN SORTING DONE. ELLEN CAN RESORT SHORTER FILES IF NEEDED LATER
    %'descriptiveFilename', '06-094327-001_to_011', ...
changeMe = changeMe + 1;
dataInfo(changeMe) = struct( ...
    'subject', 'RIHE1', ...
    'cellsFilename', '', ... % Sorted Units. the name of the NEX (NEUROEXPLORER) file that contains the final clusters to use for this dataset
    'descriptiveFilename', 'RIHE1_06-094327-001_to_011', ...
    'rethresholded', 1, ... % 0=clustered as original NEX, 1=rethresholded (with or without concatenation... this has 1 channel per file. WERE marked as INTs by Jason: 36,47,49,70,95 - I changed these all to 0 for this paper. They may all be NFS cells
    'cellChannel',       [6 13 13 13 13 19 19 20 22 22 28 28 28 29 30 30 31 31 32 32 34 36 36 36 36 41 41 41 42 42 43 45 46 46 47 47 47 49 49 51 51 52 53 53 54 54 60 60 61 62 62 67 68 68 68 68 69 70 70 71 71 74 75 78 78 86 87 87 89 90 91 95 95 96], ...
    'cellIsInterneuron', [], ...
    'cellQuality',       [4  4  3  2  1  3  3  3  4  2  4  3  2  4  3  1  4  2  4  2  3  3  4  4  3  4  3  2  4  3  3  3  4  3  3  4  2  4  3  4  3  3  4  3  4  2  5  3  4  4  3  4  5  4  3  2  4  4  2  3  2  4  4  4  3  3  4  3  4  5  4  3  4  3], ...
    ... % Seizure Info:
    'numSz', 2, ... % the number of seizures in this session
    'szStartTime', [732 19884], ...
    'szEndTime', [732+tjEndAdj_RIHE1_Sz1 19958], ...
    'szPlotStartTime', 0, ...
    'szPlotEndTime', 0, ...
    ... % Matching NS5 Info:
    'numNS5Files', 11, ... % the number of LFP files that are to be analyzed as one
    'ns5Filenames', {{'06-094327-001', '06-094327-002', '06-094327-003', '06-094327-004', ...
                      '06-094327-005', '06-094327-006', '06-094327-007', '06-094327-008', ...
                      '06-094327-009', '06-094327-010', '06-094327-011'}}, ...
    'ns5VisChannels', [1:96], ... % As seen on Blackrock Central. Leave blank to use all
    'ns5SeqChannels', [], ... % Sequential channels, 
    'ns5TrigVisChannel', 129, ... % The trigger channel, as seen on Blackrock Central
    'ns5TrigSeqChannel', 97, ...  % The sequential number for the trigger channel
    ... % Matching ECoG Info:
    'numEDFFiles', 3, ... % set to 0 if no corresponding ECoG files, or bad triggers, or BrainGate data. EcoG filenames: one entry for each NS5 that it matches
    'edfFilenames', {{'', '', ''}}, ... 
    'edfClosestVisChannels', {{''}}, ...
    'edfChannels', [], ... 
    'edfTrigVisChannel', 'Event', ...
    'edfTrigSeqChannel', 1 ...
);

%% id = 7 MG29 SZ22 all channels - ELLEN HAS SORTED ONE OR TWO CHANNELS
changeMe = changeMe + 1;
dataInfo(changeMe) = struct( ...
    'subject', 'MG29', ...
    'cellsFilename', '', ... % 20090212-111822-101-04.nex, Include suffix. Sorted Units. the name of the NEX (NEUROEXPLORER) file that contains the final clusters to use for this dataset
    'descriptiveFilename', 'MG29_seizure22_20090212-111822-100_101_102', ... % This is just a shortform for the the cellsFilename when saving heavily processed data
    'rethresholded', 1, ... % 0=clustered as original NEX, 1=rethresholded (with or without concatenation... this has 1 channel per file
    'cellChannel',       [92], ...
    'cellIsInterneuron', [1], ...
    'cellQuality',       [5], ...  
    ... % Seizure Info:
    'numSz', 1, ...
    'szStartTime', (1860 + 2400), ... % NOT EXACT
    'szEndTime', 0, ...
    'szPlotStartTime', (1840 + 2400), ...
    'szPlotEndTime', (1960 + 2400), ...
    ... % Matching NS5 Info:
    'numNS5Files', 3, ... % the number of LFP files that are to be analyzed as one
    'ns5Filenames', {{'20090212-111822-100', '20090212-111822-101', '20090212-111822-102'}}, ... % a cell array
    'ns5VisChannels', [1:96], ... % As seen on Blackrock Central. Leave blank to use all
    'ns5SeqChannels', [], ... % Sequential channels, 
    'ns5TrigVisChannel', 129, ... % The trigger channel, as seen on Blackrock Central
    'ns5TrigSeqChannel', 97, ...  % The sequenctial number for the trigger channel
    ... % Matching ECoG Info:
    'numEDFFiles', 3, ... % set to 0 if no corresponding ECoG files, or bad triggers, or BrainGate data. EcoG filenames: one entry for each NS5 that it matches
    'edfFilenames', {{'', '', ''}}, ... 
    'edfClosestVisChannels', {{''}}, ...
    'edfChannels', [], ...
    'edfTrigVisChannel', 'Event', ...
    'edfTrigSeqChannel', 1 ...
);

