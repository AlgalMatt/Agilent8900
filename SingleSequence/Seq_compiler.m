%This script requires seqauto.m, seqman1.m and seqman2. Be sure that they
%are in the current directory.

%This script takes all the data from a single Agilent 8900 sequence and
%processes it using the standard carbonate processing method. For the
%script to work you must have the entire output folder produced by the
%Agilent 8900. It can be located anywhere on your device. When running the
%script you are asked to select the BatchLog.csv file that is always 
%produced by the Agilent 8900 in the sequence output folder. This helps the
%rest of the script find said folder. 

%Automatic finds the blank and STGTE by scanning the sample names.

%Manual lets you choose where the blanks and STGTE are located. Choose this
%option only if your blank names do not contain "blk" (not
% case sensitive) somewhere in their name or the STGTE standards do not
% contain "stgte" not case sensitive) somewhere in their name. Can be used
% to correct to standards that are not STGTE.

%The data are exported to the sequence folder as an excel file. 
%After running this script you can run the Seq_plot script to quickly make
%some plots.


%Notes:
%If data is NaN, it is probably because its signal intensity is below the 
%limit of detection (here defined as 3*blank). 
%
%Blank correction and standard (STGTE) correction are both performed by
%2-point linear regression. If the sample isn't bracketted by a standard or
%blank then it will use a one-point correction. 

clear all
[~,path,~] = uigetfile('.csv', 'Select the BatchLog file in the sequence folder you want processed');
[automan,~] = listdlg('ListString',["Automatic (recommended)", "Manual"], ...
    'Name', 'Auto or manual',  'SelectionMode','single', 'ListSize', [400 400]);

if automan==1
    [CPS_t, raw_SD_t, ~, ~, bcorr_t, bcorr_SD_t, Ca43r_t, Ca43r_SD_t, STGcorr_t, STGcorr_SD_t]=seqauto(path);
else
    [CPS_t, raw_SD_t, ~, ~]=seqman1(path);
    
    [blkrows,~] = listdlg('ListString',CPS_t.Sample, 'Name', 'Select the blanks', 'ListSize', [400 400]);
    [STGTEindx,~] = listdlg('ListString',CPS_t.Sample, 'Name', 'Select the STGTE', 'ListSize', [400 400]);
    
    [bcorr_t, bcorr_SD_t, Ca43r_t, Ca43r_SD_t, STGcorr_t, STGcorr_SD_t]=seqman2(CPS_t,...
        raw_SD_t, blkrows, STGTEindx);    
end

DatTab{1}=CPS_t; ErrTab{1}=raw_SD_t;
DatTab{2}=bcorr_t; ErrTab{2}=bcorr_SD_t;
DatTab{3}=Ca43r_t; ErrTab{3}=Ca43r_SD_t;
DatTab{4}=STGcorr_t; ErrTab{4}=STGcorr_SD_t;

Elements=CPS_t.Properties.VariableNames(4:end);

save('Seq_Ts.mat', 'DatTab', 'ErrTab', 'Elements');

slashIdx = strfind(path, '\');
runname=path(slashIdx(end-1)+1:slashIdx(end)-3);
writetable(CPS_t,[path, runname, '.xlsx'],'FileType','spreadsheet','Sheet','Raw CPS')
writetable(raw_SD_t,[path, runname, '.xlsx'],'FileType','spreadsheet','Sheet','Raw CPS 1sd')
writetable(bcorr_t,[path, runname, '.xlsx'],'FileType','spreadsheet','Sheet','Blank corrected')
writetable(bcorr_SD_t,[path, runname, '.xlsx'],'FileType','spreadsheet','Sheet','Blank corr 1sd')
writetable(Ca43r_t,[path, runname, '.xlsx'],'FileType','spreadsheet','Sheet','Ca43 ratios')
writetable(Ca43r_SD_t,[path, runname, '.xlsx'],'FileType','spreadsheet','Sheet','Ca43 1sd')
writetable(STGcorr_t,[path, runname, '.xlsx'],'FileType','spreadsheet','Sheet','STGTE corrected')
writetable(STGcorr_SD_t,[path, runname, '.xlsx'],'FileType','spreadsheet','Sheet','STG corr 1sd')
