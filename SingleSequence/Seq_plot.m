%This script requires a Seq_Ts.mat file (produced by running the
%Seq_autocompiler.m or Seq_mancompiler.m) and a Consistency.mat file, which
%you can get from me or by running the Agilent_process.m.
%This script produces quick plots of the data from a Agilent 8900 sequence.

%Simple press Run and click the options.

clear all
load 'Seq_Ts.mat'
load 'Consistency.mat'

[Pindx,~] = listdlg('ListString',["Raw CPS", "Blank corrected", "Me/Ca ratio",...
    "STG-corrected", "Consistency standards", "Blanks as fraction of STGTE"],...
    'Name', 'Select type of data you want to plot', 'SelectionMode','single', 'ListSize', [400 400]);

if Pindx < 5
 [Sindx,~] = listdlg('ListString',DatTab{Pindx}.Sample,...
    'Name', 'Select samples you want to plot', 'ListSize', [400 400]);
end
[Eindx,~] = listdlg('ListString',Elements,...
    'Name', 'Select elements you want to plot', 'ListSize', [400 400]);


%% Figures
if any(Pindx==1)
for i=1:numel(Eindx) 
    %raw CPS
f1=figure('units','normalized','outerposition',[0 0 1 1]);
errorbar(datenum(DatTab{1}{Sindx, 'Time'}), DatTab{1}{Sindx, Elements{Eindx(i)}}, 2*ErrTab{1}{Sindx, Elements{Eindx(i)}}, 'o-r')
datetick('x',13)
title([Elements{Eindx(i)} ' raw counts ' char(177) '2\sigma'])
xlabel('Time')
ylabel('Counts per second')
 text(datenum(DatTab{1}{Sindx, 'Time'}),DatTab{1}{Sindx, Elements{Eindx(i)}},DatTab{1}.Sample(Sindx),'VerticalAlignment','bottom','HorizontalAlignment','right')
grid on
end
end

if any(Pindx==2)
for i=1:numel(Eindx) 
    %bcorr
f2=figure('units','normalized','outerposition',[0 0 1 1]);
errorbar(datenum(DatTab{2}{Sindx, 'Time'}), DatTab{2}{Sindx, Elements{Eindx(i)}}, 2*ErrTab{2}{Sindx, Elements{Eindx(i)}}, 'o-r')
datetick('x',13)
title([Elements{Eindx(i)} ' blank-corrected ' char(177) '2\sigma'])
xlabel('Time')
ylabel('Counts per second')
 text(datenum(DatTab{2}{Sindx, 'Time'}),DatTab{2}{Sindx, Elements{Eindx(i)}},DatTab{2}.Sample(Sindx),'VerticalAlignment','bottom','HorizontalAlignment','right')
grid on
end
end

if any(Pindx==3)
for i=1:numel(Eindx) 
    %Me/Ca
f3=figure('units','normalized','outerposition',[0 0 1 1]);
errorbar(datenum(DatTab{3}{Sindx, 'Time'}), DatTab{3}{Sindx, Elements{Eindx(i)}}, 2*ErrTab{3}{Sindx, Elements{Eindx(i)}}, 'o-r')
datetick('x',13)
title([Elements{Eindx(i)}, '/Ca43 ' char(177) '2\sigma'])
xlabel('Time')
ylabel([Elements{Eindx(i)}, '/Ca43 (CPS/CPS)'])
 text(datenum(DatTab{3}{Sindx, 'Time'}),DatTab{3}{Sindx, Elements{Eindx(i)}},DatTab{3}.Sample(Sindx),'VerticalAlignment','bottom','HorizontalAlignment','right')
grid on
end
end

if any(Pindx==4)
   %Elements=Elements(~contains(Elements, 'Ca'));
for i=1:numel(Eindx) 
    %Me/Ca
f4=figure('units','normalized','outerposition',[0 0 1 1]);
errorbar(datenum(DatTab{4}{Sindx, 'Time'}), DatTab{4}{Sindx, Elements{Eindx(i)}}, 2*ErrTab{4}{Sindx, Elements{Eindx(i)}}, 'o-r')
datetick('x',13)
title([Elements{Eindx(i)}, '/Ca43 ' char(177) '2\sigma'])
xlabel('Time')
ylabel([Elements{Eindx(i)}, '/Ca43 (\mumol/mol)'])
 text(datenum(DatTab{4}{Sindx, 'Time'}),DatTab{4}{Sindx, Elements{Eindx(i)}},DatTab{4}.Sample(Sindx),'VerticalAlignment','bottom','HorizontalAlignment','right')
grid on
end
end


if any(Pindx==5)
CS1idx=contains(lower(DatTab{1}.Sample), 'cs1'); 
CS2idx=contains(lower(DatTab{1}.Sample), 'cs2'); 
CS3idx=contains(lower(DatTab{1}.Sample), 'cs3');
for i=1:numel(Eindx) 
    %Me/Ca
f5=figure('units','normalized','outerposition',[0 0 1 1]);
p1=plot(CS1_T.Time, CS1_T{:, Elements{Eindx(i)}}, 'ob');
yline(nanmean(CS1_T{:, Elements{Eindx(i)}}), 'b');
yline(nanmean(CS1_T{:, Elements{Eindx(i)}})+nanstd(CS1_T{:, Elements{Eindx(i)}})*2, '-.b');
yline(nanmean(CS1_T{:, Elements{Eindx(i)}})-nanstd(CS1_T{:, Elements{Eindx(i)}})*2, '-.b');
hold on
p2=plot(DatTab{4}{CS1idx, 'Time'}, DatTab{4}{CS1idx, Elements{Eindx(i)}}, 'or', 'MarkerFaceColor', 'r');
title(['CS1 ',Elements{Eindx(i)}, '/Ca43 ' char(177) '2\sigma'])
xlabel('Time')
ylabel([Elements{Eindx(i)}, '/Ca43 (\mumol/mol)'])
 %text(datenum(DatTab{4}{Sindx, 'Time'}),DatTab{4}{Sindx, Elements{Eindx(i)}},DatTab{4}.Sample(Sindx),'VerticalAlignment','bottom','HorizontalAlignment','right')
grid on
legend([p1, p2], 'Previous data', 'This run')
hold off

f6=figure('units','normalized','outerposition',[0 0 1 1]);
p1=plot(CS2_T.Time, CS2_T{:, Elements{Eindx(i)}}, 'ob');
yline(nanmean(CS2_T{:, Elements{Eindx(i)}}), 'b');
yline(nanmean(CS2_T{:, Elements{Eindx(i)}})+nanstd(CS2_T{:, Elements{Eindx(i)}})*2, '-.b');
yline(nanmean(CS2_T{:, Elements{Eindx(i)}})-nanstd(CS2_T{:, Elements{Eindx(i)}})*2, '-.b');
hold on
p2=plot(DatTab{4}{CS2idx, 'Time'}, DatTab{4}{CS2idx, Elements{Eindx(i)}}, 'or', 'MarkerFaceColor', 'r');
title(['CS2 ', Elements{Eindx(i)}, '/Ca43 ' char(177) '2\sigma'])
xlabel('Time')
ylabel([Elements{Eindx(i)}, '/Ca43 (\mumol/mol)'])
 %text(datenum(DatTab{4}{Sindx, 'Time'}),DatTab{4}{Sindx, Elements{Eindx(i)}},DatTab{4}.Sample(Sindx),'VerticalAlignment','bottom','HorizontalAlignment','right')
grid on
legend([p1, p2], 'Previous data', 'This run')
hold off

f7=figure('units','normalized','outerposition',[0 0 1 1]);
p1=plot(CS3_T.Time, CS3_T{:, Elements{Eindx(i)}}, 'ob');
yline(nanmean(CS3_T{:, Elements{Eindx(i)}}), 'b');
yline(nanmean(CS3_T{:, Elements{Eindx(i)}})+nanstd(CS3_T{:, Elements{Eindx(i)}})*2, '-.b');
yline(nanmean(CS3_T{:, Elements{Eindx(i)}})-nanstd(CS3_T{:, Elements{Eindx(i)}})*2, '-.b');
hold on
p2=plot(DatTab{4}{CS3idx, 'Time'}, DatTab{4}{CS3idx, Elements{Eindx(i)}}, 'or', 'MarkerFaceColor', 'r');
title(['CS3 ', Elements{Eindx(i)}, '/Ca43 ' char(177) '2\sigma'])
xlabel('Time')
ylabel([Elements{Eindx(i)}, '/Ca43 (\mumol/mol)'])
grid on
legend([p1, p2], 'Previous data', 'This run')
hold off
end
end


if any(Pindx==6)
    
    %find locations of STGTE
STGTErows=find(contains(lower(DatTab{1}.Sample),'stgte'));

%Discount any 0.5STGTE
if ~all(contains(lower(DatTab{1}.Sample(STGTErows)), '0.5stgte'))
    STGTErows(contains(lower(DatTab{1}.Sample(STGTErows)), '0.5stgte'))=[];
end
   
for i=1:numel(Eindx) 
    %blanks calibrated to STGTE
y=DatTab{1}{contains(lower(DatTab{1}.Sample),'blk'), Elements{Eindx(i)}}./nanmean(DatTab{1}{STGTErows,Elements{Eindx(i)}},1) ;
x=DatTab{1}{contains(lower(DatTab{1}.Sample),'blk'), 'Time'};
f8=figure('units','normalized','outerposition',[0 0 1 1]);
p1=plot(bmean_cali_T.Time,bmean_cali_T{:,Elements{Eindx(i)}},'ob');
hold on
yline(nanmean(bmean_cali_T{:,Elements{Eindx(i)}}), 'b');
yline(nanmean(bmean_cali_T{:,Elements{Eindx(i)}})+nanstd(bmean_cali_T{:,Elements{Eindx(i)}})*2, '-.b');
yline(nanmean(bmean_cali_T{:,Elements{Eindx(i)}})-nanstd(bmean_cali_T{:,Elements{Eindx(i)}})*2, '-.b');
p2=plot(x, y, 'or', 'MarkerFaceColor', 'r');
title(['STGTE-calibrated blanks ' Elements{Eindx(i)}  ' ' char(177) '2\sigma'])
xlabel('Time')
ylabel('Blank CPS as fraction of STGTE CPS')
grid on
legend([p1, p2], 'Previous data', 'This run')
hold off
end
end



