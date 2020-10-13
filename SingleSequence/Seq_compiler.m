%This script requires seqauto.m, seqman1.m and seqman2. For the figures
%you'll need Consistency.mat. Be sure that they are in the current 
%directory.

%This script takes all the data from a single Agilent 8900 sequence and
%processes it using the normal carbonate processing method. For the
%script to work you must have the entire output folder produced by the
%Agilent 8900 (not just the summary file). It can be located anywhere on 
%your device. When running the script you are asked to select the 
%BatchLog.csv file that is always produced by the Agilent 8900 in the 
%sequence output folder. This helps the rest of the script find said 
%folder. 

%You are given two options:
%Automatic finds the blank and STGTE by scanning the sample names. If the
%sample names have STGTE in any form (i.e. stgte, 1stgte, STGTE etc.), 
%they will be identified as the STGTE standards. Same goes for blanks 
%(labelled as BLK in any form).

%Manual lets you choose where the blanks and STGTE are located. Can be used
%to correct to standards that are not STGTE or omit certain stnds or blks.

%The data are exported to the sequence folder as an excel file. 

%There is a further option to make quick plots. This will produce a
%selection of plots of the data as pdf files exported to the sequence
%folder. Use them as a quick check of the data. 

%The plots include...
%STGcorr: the final data from all the samples and consistency standards.

%CS1, CS2, CS3, 8301f: Any consistency standards run plotted with a
%database of past values of the consistency standards. Good to check.

%Blanks: The blanks from the run plotted with mean blanks from previous
%runs normalised to STGTE. 

%After running this script you can also run the Seq_plot script to quickly 
%make some plots.

%Notes:  
%Blank correction and standard (STGTE) correction are both performed by
%2-point linear regression between the closest two STGTE or Blks. If the 
%sample isn't bracketted by a standard or blank then it will use a 
%one-point correction. 

clear all
[~,path,~] = uigetfile('.csv', 'Select the BatchLog file in the sequence folder you want processed');
[automan,~] = listdlg('ListString',["Automatic (recommended)", "Manual"], ...
    'Name', 'Automatic or manual compilin',  'SelectionMode','single', 'ListSize', [400 400]);

[makefigs,~] = listdlg('ListString',["Yes", "No"], ...
    'Name', 'Make quick figures?',  'SelectionMode','single', 'ListSize', [400 400]);

disp('Compiling...')
if automan==1
    [CPS_t, raw_SD_t, ~, ~, bcorr_t, bcorr_SD_t, Ca43r_t, Ca43r_SD_t, STGcorr_t, STGcorr_SD_t, STGTErows]=seqauto(path);
else
    [CPS_t, raw_SD_t, ~, ~]=seqman1(path);
    
    [blkrows,~] = listdlg('ListString',CPS_t.Sample, 'Name', 'Select the blanks', 'ListSize', [400 400]);
    [STGTErows,~] = listdlg('ListString',CPS_t.Sample, 'Name', 'Select the STGTE', 'ListSize', [400 400]);
    
    [bcorr_t, bcorr_SD_t, Ca43r_t, Ca43r_SD_t, STGcorr_t, STGcorr_SD_t]=seqman2(CPS_t,...
        raw_SD_t, blkrows, STGTErows);    
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
disp('Compilation complete. Exported output to:')
disp([path, runname, '.xlsx'])





%% Figures 
if makefigs ==1
disp('Carefully sketching figures...')
%Stg-corrected all samples
STGTE_elements={'Li7', 'B11', 'Na23', 'Mg24', 'Mg25', 'Al27', 'Mn55', ...
    'Sr88', 'Cd111', 'Ba138', 'Nd146', 'U238'};
panfig1=figure('units','centimeters');
pos = get(panfig1,'Position');
set(panfig1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[29.7, 42])
for i=1:numel(STGTE_elements) 
 subplot(ceil(numel(STGTE_elements)/2), 2, i)   
errorbar(datenum(DatTab{4}{:, 'Time'}), DatTab{4}{:, STGTE_elements{i}},...
    1*ErrTab{4}{:, STGTE_elements{i}}, 'o-r', 'MarkerSize', 1, 'CapSize', 2)
datetick('x',13)
ax=gca;
ax.XAxis.FontSize=2;
ax.YAxis.FontSize=2;
title([STGTE_elements{i}, '/Ca43 ' char(177) '1\sigma'], 'FontSize', 8)
xlabel('Time', 'FontSize', 6)
ylabel([STGTE_elements{i}, '/Ca43 (\mumol/mol)'], 'FontSize', 8)
 text(datenum(DatTab{4}{:, 'Time'}),DatTab{4}{:, STGTE_elements{i}},...
     DatTab{4}.Sample,'VerticalAlignment','bottom','HorizontalAlignment','right', 'FontSize',8)
 grid on
hold off
sgtitle('STGTE-corrected Element/Ca43 ratios')
end
print(panfig1, [path,'STGcorr'] ,'-dpdf', '-fillpage')

%CS
load 'Consistency.mat'
CS1idx=contains(lower(DatTab{1}.Sample), 'cs1');
CS2idx=contains(lower(DatTab{1}.Sample), 'cs2');
CS3idx=contains(lower(DatTab{1}.Sample), 'cs3');
s8301fidx=contains(lower(DatTab{1}.Sample), '8301f');
%CS1
panfig2=figure('units','centimeters');
pos = get(panfig2,'Position');
set(panfig2,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[29.7, 42])
for i=1:numel(STGTE_elements)
    subplot(ceil(numel(STGTE_elements)/2), 2, i)
    p1=plot(CS1_T.Time, CS1_T{:, STGTE_elements{i}}, 'ob', 'MarkerSize', 4);
    yline(nanmean(CS1_T{:, STGTE_elements{i}}), 'b');
    yline(nanmean(CS1_T{:, STGTE_elements{i}})+nanstd(CS1_T{:, STGTE_elements{i}}), '-.b');
    yline(nanmean(CS1_T{:, STGTE_elements{i}})-nanstd(CS1_T{:, STGTE_elements{i}}), '-.b');
    hold on
    p2=plot(DatTab{4}{CS1idx, 'Time'}, DatTab{4}{CS1idx, STGTE_elements{i}}, 'or', 'MarkerFaceColor', 'r', 'MarkerSize', 4);
    title(['CS1 ',STGTE_elements{i}, '/Ca43 ' char(177) '1\sigma'])
    xlabel('Time')
    ylabel([STGTE_elements{i}, '/Ca43 (\mumol/mol)'])
    grid on
  legend([p1, p2], 'Previous data', 'This run', 'FontSize', 8, 'Position',[0.58 0.915 0.02 0.05])
    legend('boxoff')
    sgtitle('CS1')
    hold off
end
print(panfig2, [path,'CS1'] ,'-dpdf', '-fillpage')      

%CS2
panfig3=figure('units','centimeters');
pos = get(panfig3,'Position');
set(panfig3,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[29.7, 42])
for i=1:numel(STGTE_elements)
    subplot(ceil(numel(STGTE_elements)/2), 2, i)
    p1=plot(CS2_T.Time, CS2_T{:, STGTE_elements{i}}, 'ob', 'MarkerSize', 4);
    yline(nanmean(CS2_T{:, STGTE_elements{i}}), 'b');
    yline(nanmean(CS2_T{:, STGTE_elements{i}})+nanstd(CS2_T{:, STGTE_elements{i}}), '-.b');
    yline(nanmean(CS2_T{:, STGTE_elements{i}})-nanstd(CS2_T{:, STGTE_elements{i}}), '-.b');
    hold on
    p2=plot(DatTab{4}{CS2idx, 'Time'}, DatTab{4}{CS2idx, STGTE_elements{i}}, 'or', 'MarkerFaceColor', 'r', 'MarkerSize', 4);
    title(['CS1 ',STGTE_elements{i}, '/Ca43 ' char(177) '1\sigma'])
    xlabel('Time')
    ylabel([STGTE_elements{i}, '/Ca43 (\mumol/mol)'])
    grid on
      legend([p1, p2], 'Previous data', 'This run', 'FontSize', 8, 'Position',[0.58 0.915 0.02 0.05])
    legend('boxoff')
    sgtitle('CS2')
    hold off
end
print(panfig3, [path, 'CS2'] ,'-dpdf', '-fillpage')   


%CS3
panfig4=figure('units','centimeters');
pos = get(panfig4,'Position');
set(panfig4,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[29.7, 42])
for i=1:numel(STGTE_elements)
    subplot(ceil(numel(STGTE_elements)/2), 2, i)
    p1=plot(CS3_T.Time, CS3_T{:, STGTE_elements{i}}, 'ob', 'MarkerSize', 4);
    yline(nanmean(CS3_T{:, STGTE_elements{i}}), 'b');
    yline(nanmean(CS3_T{:, STGTE_elements{i}})+nanstd(CS3_T{:, STGTE_elements{i}}), '-.b');
    yline(nanmean(CS3_T{:, STGTE_elements{i}})-nanstd(CS3_T{:, STGTE_elements{i}}), '-.b');
    hold on
    p2=plot(DatTab{4}{CS3idx, 'Time'}, DatTab{4}{CS3idx, STGTE_elements{i}}, 'or', 'MarkerFaceColor', 'r', 'MarkerSize', 4);
    title(['CS1 ',STGTE_elements{i}, '/Ca43 ' char(177) '1\sigma'])
    xlabel('Time')
    ylabel([STGTE_elements{i}, '/Ca43 (\mumol/mol)'])
    grid on
    legend([p1, p2], 'Previous data', 'This run', 'FontSize', 8, 'Position',[0.58 0.915 0.02 0.05])
    legend('boxoff')
    sgtitle('CS3')
    hold off
end
print(panfig4, [path, 'CS3'] ,'-dpdf', '-fillpage')  

%8301f
panfig5=figure('units','centimeters');
pos = get(panfig5,'Position');
set(panfig5,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[29.7, 42])
for i=1:numel(STGTE_elements)
    subplot(ceil(numel(STGTE_elements)/2), 2, i)
    p1=plot(std8301f_T.Time, std8301f_T{:, STGTE_elements{i}}, 'ob', 'MarkerSize', 4);
    yline(nanmean(std8301f_T{:, STGTE_elements{i}}), 'b');
    yline(nanmean(std8301f_T{:, STGTE_elements{i}})+nanstd(std8301f_T{:, STGTE_elements{i}}), '-.b');
    yline(nanmean(std8301f_T{:, STGTE_elements{i}})-nanstd(std8301f_T{:, STGTE_elements{i}}), '-.b');
    hold on
    p2=plot(DatTab{4}{s8301fidx, 'Time'}, DatTab{4}{s8301fidx, STGTE_elements{i}}, 'or', 'MarkerFaceColor', 'r', 'MarkerSize', 4);
    title(['CS1 ',STGTE_elements{i}, '/Ca43 ' char(177) '1\sigma'])
    xlabel('Time')
    ylabel([STGTE_elements{i}, '/Ca43 (\mumol/mol)'])
    grid on
    legend([p1, p2], 'Previous data', 'This run', 'FontSize', 8, 'Position',[0.58 0.915 0.02 0.05])
    legend('boxoff')
    sgtitle('8301f')
    hold off
   
end
print(panfig5, [path, '8301f'] ,'-dpdf', '-fillpage')


%Blanks
panfig6=figure('units','centimeters');
pos = get(panfig6,'Position');
set(panfig6,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[29.7, 42])
for i=1:numel(STGTE_elements)
    subplot(ceil(numel(STGTE_elements)/2), 2, i)
    y=DatTab{1}{contains(lower(DatTab{1}.Sample),'blk'),...
        STGTE_elements{i}}./nanmean(DatTab{1}{STGTErows,STGTE_elements{i}},1) ;
    x=DatTab{1}{contains(lower(DatTab{1}.Sample),'blk'), 'Time'};
    p1=plot(bmean_cali_T.Time,bmean_cali_T{:,STGTE_elements{i}},'ob', 'MarkerSize', 4);
    hold on
    yline(nanmean(bmean_cali_T{:,STGTE_elements{i}}), 'b');
    yline(nanmean(bmean_cali_T{:,STGTE_elements{i}})+nanstd(bmean_cali_T{:,STGTE_elements{i}}), '-.b');
    yline(nanmean(bmean_cali_T{:,STGTE_elements{i}})-nanstd(bmean_cali_T{:,STGTE_elements{i}}), '-.b');
    p2=plot(x, y, 'or', 'MarkerFaceColor', 'r', 'MarkerSize', 4);
    title(['STGTE-calibrated blanks ' STGTE_elements{i}  ' ' char(177) '1\sigma'])
    xlabel('Time')
    ylabel('Blank CPS as fraction of STGTE CPS')
    grid on
    legend([p1, p2], 'Previous data', 'This run', 'FontSize', 8, 'Position',[0.58 0.915 0.02 0.05])
    legend('boxoff')
    sgtitle('Blanks')
    hold off
end

print(panfig6, [path, 'Blanks'] ,'-dpdf', '-fillpage')

disp('Figures done.')
disp('Exported output to...')
disp([path, runname, '.xlsx'])
end

disp('Woohoo! All done.')

