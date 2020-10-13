function [CPS_t, raw_SD_t, raw_N_t, raw_intTime_t, bcorr_t, bcorr_SD_t, Ca43r_t, Ca43r_SD_t, STGcorr_t, STGcorr_SD_t, STGTErows] = seqauto(path)

Batchlogloc=fullfile(path,'BatchLog.csv');

optsbatch = detectImportOptions(Batchlogloc);
optsbatch = setvaropts(optsbatch,'Acq_Date_Time','InputFormat','MM/dd/uuuu hh:mm:ss aa');
batchlog=readtable(Batchlogloc,optsbatch);
Samplelist=batchlog.SampleName(ismember(batchlog.AcquisitionResult,'Pass'));
slashIdx = strfind(path, '\');
runname=path(slashIdx(end-1)+1:slashIdx(end)-1);

filelist=cell2mat(batchlog.FileName(ismember(batchlog.AcquisitionResult,'Pass') & ~ismember(batchlog.FileName,'-')));
subFolders=string(filelist(:,end-8:end));

for i = 1:size(subFolders,1)
    
    file=dir(fullfile(path,subFolders(i),'*1.csv'));
    if isempty(file)
        continue
    end
    Fileloc=fullfile(file.folder,file.name);
    
    %Sets import options for CSV file. Which rows to import and where the
    %variable names are located.
    opts = detectImportOptions(Fileloc, 'Delimiter', ',', 'NumHeaderLines', 7);
    opts.Delimiter=',';
    opts.VariableNamesLine = 8;
    opts.DataLines = [9 inf];
    opts.RowNamesColumn=1;
    
    %import data
    rawdat=readtable(Fileloc, opts);
    %get rid of excess info at end of file
    printinfo=find(contains(rawdat.Properties.RowNames,'Print'));
    rawdat([printinfo],:)=[];
    %extract data from individual run output file
    RSDs=rawdat.RSD___; %RSDs
    RSDs(RSDs==0)=NaN;
    SDs=rawdat.SD;
    SDs(SDs==0)=NaN;
    
    Samplename=string(Samplelist(i));
    Time=datenum(batchlog.Acq_Date_Time(i));
    Time=datetime(Time, 'Format', 'yyyy-MM-dd HH:mm:ss', 'convertFrom', 'datenum');
    
    SampleData=[{runname,Time,Samplename}, num2cell(rawdat.CPS)'];
    SampleSDs=[{runname,Time,Samplename}, num2cell(SDs)'];
    SampleN=[{runname,Time,Samplename}, num2cell(rawdat.n)'];
    SampleintTime=[{runname,Time,Samplename}, num2cell(rawdat.Time_Sec_)'];
    
    Elements=strcat(rawdat.Element,(""+rawdat.Mass.')');
    
    Titles=['RunName','Time','Sample',Elements'];
    
    SampleData_t=cell2table(SampleData,'VariableNames', Titles);
    SampleSDs_t=cell2table(SampleSDs,'VariableNames', Titles);
    SampleN_t=cell2table(SampleN,'VariableNames', Titles);
    SampleintTime_t=cell2table(SampleintTime,'VariableNames', Titles);
    
    %put all of the sample tables together to make a run table
    if i == 1
        run_CPS_t=SampleData_t;
        run_SDs_t=SampleSDs_t;
        run_N_t=SampleN_t;
        run_intTime_t=SampleintTime_t;
    else
        
        %returns the elements that are present in the sample file but
        %not in the rest of the run (until this sample). All samples
        %should be measured identically through the run anyway so this
        %code may be redundant.
        a=setdiff(SampleData_t.Properties.VariableNames,run_CPS_t.Properties.VariableNames);
        %returns the elements that are present in the run but
        %not in this sample.
        b=setdiff(run_CPS_t.Properties.VariableNames,SampleData_t.Properties.VariableNames);
        
        %puts nan in places where there are differences in elements between
        %sample and run.
        t2_data = [SampleData_t array2table(nan(height(SampleData_t), numel(b)), 'VariableNames', b)];
        T2_data= [run_CPS_t array2table(nan(height(run_CPS_t), numel(a)), 'VariableNames', a)];
        
        t2_SDs=[SampleSDs_t array2table(nan(height(SampleSDs_t), numel(b)), 'VariableNames', b)];
        T2_SDs=[run_SDs_t array2table(nan(height(run_SDs_t), numel(a)), 'VariableNames', a)];
        
        t2_N=[SampleN_t array2table(nan(height(SampleN_t), numel(b)), 'VariableNames', b)];
        T2_N=[run_N_t array2table(nan(height(run_N_t), numel(a)), 'VariableNames', a)];
        
        t2_intTime=[SampleintTime_t array2table(nan(height(SampleintTime_t), numel(b)), 'VariableNames', b)];
        T2_intTime=[run_intTime_t array2table(nan(height(run_intTime_t), numel(a)), 'VariableNames', a)];
        
        %Puts sample into run table
        run_CPS_t=[T2_data;t2_data];
        run_SDs_t=[T2_SDs;t2_SDs];
        run_N_t=[T2_N;t2_N];
        run_intTime_t=[T2_intTime;t2_intTime];
        
    end
    
end %%%%%%%%%%%%%end of cycling through samples in a run

CPS_t=sortrows(run_CPS_t,{'Time'},{'ascend'});
raw_SD_t=sortrows(run_SDs_t,{'Time'},{'ascend'});
raw_N_t=sortrows(run_N_t,{'Time'},{'ascend'});
raw_intTime_t=sortrows(run_intTime_t,{'Time'},{'ascend'});



%% Processing


Elements=CPS_t.Properties.VariableNames(4:end); %Pulls full list of elements

%theoretical error vs real error
%N_array=raw_N_t{:, Elements}; %number of cycles

%Published STGTE values
STGTE_values=[20.711, 86.1, 8.09, 5.382, 5.382, 42.382, 58.634, 1.488,...
    0.14, 1.277, 5.629, 159.567];
STGTE_elements={'Li7', 'B11', 'Na23', 'Mg24', 'Mg25', 'Al27', 'Mn55', ...
    'Sr88', 'Cd111', 'Ba138', 'Nd146', 'U238'};

elapse=CPS_t.Time-CPS_t.Time(1); %Time since beginning of run
elapse(elapse>duration(60,00,00))=duration(00,00,00); %max limit on elapse. Stops duplicate analyses.
info=[table2cell(CPS_t(:,1:3)), num2cell(elapse)];
realSD_t=raw_SD_t(:,Elements);
%run_N=N_array;

%automatic blanks
blkrows=find(contains(lower(CPS_t.Sample), 'blk')); %blank locations
%% blank corrections

blks=CPS_t{blkrows,Elements};  %blank data only (needs to be array for processing below)
blks_SD=realSD_t{blkrows, Elements};
%find outliers in blanks. Replace them with nan.
bol=blks>nanmean(blks,1)+2.5*nanstd(blks,1) | blks<nanmean(blks,1)-2.5*nanstd(blks,1);
blkrows_ol=blkrows;
while ~isempty(find(bol, 1))
    blks(bol)=NaN;
    blks_SD(bol)=nan;
    blkrows_ol(any(bol, 2))=nan;
    bol=blks>nanmean(blks,1)+2.5*nanstd(blks,1) | blks<nanmean(blks,1)-2.5*nanstd(blks,1);
end
blkrows_ol=blkrows_ol(~isnan(blkrows_ol));
bmean=nanmean(blks,1);
%error propagation
%bmean_SD=(nansum(blks_SD.^2,1).^0.5)./nansum(~isnan(blks_SD), 1);
bmean_t=cell2table([table2cell(CPS_t(blkrows(1),1:3)), num2cell(bmean)],...
    'VariableNames', ['RunName','Time','Sample', Elements]);
%pre-allocate tables for computational efficiency
bcorr=zeros(size(CPS_t{:,Elements}));
bcorr_SD=bcorr;

secs=seconds(elapse);
for j=1:numel(secs)
    [val, idx]=sort(blkrows_ol-j, 'ComparisonMethod', 'abs');
    %If there are not two blks bracketing the sample then just use the
    %closest single BLK as a point correction.    
    bcorr(j,:)=CPS_t{j,Elements}-CPS_t{blkrows_ol(idx(1)),Elements};
    bcorr_SD(j,:)=(realSD_t{j,Elements}.^2+realSD_t{blkrows_ol(idx(1)),Elements}.^2).^0.5;
    %Otherwise use linear regression
    if ~all(val>=0) && ~all(val<=0)
        uidx=blkrows_ol(idx(find(val>0, 1)));
        lidx=blkrows_ol(idx(find(val<=0, 1)));
        c=CPS_t{lidx, Elements};
        m=(CPS_t{uidx, Elements}-c)./(secs(uidx)-secs(lidx));
        bcorr(j,:)=CPS_t{j,Elements}-(m.*(secs(j)-secs(lidx))+c);
        
        %error propagation (regress)
        c=realSD_t{lidx, Elements};
        m=(realSD_t{uidx,Elements}-c)./(secs(uidx)-secs(lidx));
        bcorr_SD(j,:)=(realSD_t{j,Elements}.^2+(m.*(secs(j)-secs(lidx))+c).^2).^0.5;
    end
end

bcorr_RSD=bcorr_SD./bcorr;
bcorr_t=cell2table([info,num2cell(bcorr)],'VariableNames', ['RunName','Time','Sample','Elapse', Elements]);
bcorr_SD_t=cell2table([info,num2cell(bcorr_SD)],'VariableNames', ['RunName','Time','Sample','Elapse', Elements]);
%bcorr_RSD_t=cell2table([info,num2cell(bcorr_RSD)],'VariableNames', ['RunName','Time','Sample','Elapse', Elements]);

%% Me/43Ca ratios for the standards
Ca43_loc=find(ismember(Elements,'Ca43'));
Ca43r=bcorr./bcorr(:,Ca43_loc);
Ca43r_t=cell2table([info,num2cell(Ca43r)],'VariableNames', ['RunName','Time','Sample', 'Elapse', Elements]);
Ca43r_RSD=(bcorr_RSD.^2+bcorr_RSD(:, Ca43_loc).^2).^0.5;
Ca43r_RSD(isinf(Ca43r_RSD))=NaN;

%% Me/Ca counting stats
% Ca43r_RSE_t=cell2table([info,num2cell((Ca43r_RSD.*Ca43r./(run_N.^0.5))./Ca43r)],'VariableNames', ['RunName','Time','Sample', 'Elapse', Elements]);
Ca43r_RSD_t=cell2table([info,num2cell(Ca43r_RSD)],'VariableNames', ['RunName','Time','Sample', 'Elapse', Elements]);
Ca43r_SD_t=cell2table([info,num2cell(Ca43r_RSD.*Ca43r)],'VariableNames', ['RunName','Time','Sample', 'Elapse', Elements]);
%Ca43r_theoRSD_t=cell2table([info,num2cell(Ca43r_theoRSD)],'VariableNames', ['RunName','Time','Sample', Elements]);


%% standard calibration

%find locations of STGTE
STGTErows=find(contains(lower(Ca43r_t.Sample),'stgte'));

%Discount any 0.5STGTE
if ~all(contains(lower(Ca43r_t.Sample(STGTErows)), '0.5stgte'))
    STGTErows(contains(lower(Ca43r_t.Sample(STGTErows)), '0.5stgte'))=[];
end

Ca43r_stg=Ca43r_t{:,STGTE_elements};
Ca43r_RSDstg=Ca43r_RSD_t{:,STGTE_elements};

%remove anomalously low counts (blockages or exhaustion of sample)
%Ca43r_stg(bcorr_t{:,STGTE_elements}<bmean_t{:,STGTE_elements}*3 | Ca43r_stg<0)=nan;
%Ca43r_RSDstg(isnan(Ca43r_stg))=nan;

%pre-allocate tables for computational efficiency
STGcorr=zeros(size(Ca43r_stg));
STGcorr_RSD=STGcorr;

secs=seconds(Ca43r_t{:, 'Elapse'});
for j=1:numel(secs)
    [val, idx]=sort(STGTErows-j, 'ComparisonMethod', 'abs');
    %If there are not two STGTE bracketing the sample then just use the
    %closest single STGTE as a point correction.
    STGcorr(j,:)=Ca43r_stg(j,:)./Ca43r_stg(STGTErows(idx(1)),:).*STGTE_values;
    STGcorr_RSD(j,:)=(Ca43r_RSDstg(j,:).^2+Ca43r_RSDstg(STGTErows(idx(1)),:).^2).^0.5;
    
    %Otherwise use linear regression
    if ~all(val>=0) && ~all(val<=0)
        uidx=STGTErows(idx(find(val>0, 1)));
        lidx=STGTErows(idx(find(val<=0, 1)));
        c=Ca43r_stg(lidx, :);
        m=(Ca43r_stg(uidx, :)-c)./(secs(uidx)-secs(lidx));
        STGcorr(j,:)=Ca43r_stg(j,:)./(m.*(secs(j)-secs(lidx))+c).*STGTE_values;
        
        %error propagation (regress)
        c=Ca43r_RSDstg(lidx, :);
        m=(Ca43r_RSDstg(uidx,:)-c)./(secs(uidx)-secs(lidx));       
        STGcorr_RSD(j,:)=(Ca43r_RSDstg(j,:).^2+(m.*(secs(j)-secs(lidx))+c).^2).^0.5;
    end
end

%Remove outliers
STGcorr_RSD(STGcorr_RSD>1E3)=nan;

STGcorr_t=cell2table([info,num2cell(STGcorr)],'VariableNames', ['RunName','Time','Sample', 'Elapse', STGTE_elements]);
%STGcorr_RSD_t=cell2table([info,num2cell(STGcorr_RSD)],'VariableNames', ['RunName','Time','Sample','Elapse', STGTE_elements]);
STGcorr_SD_t=cell2table([info,num2cell(STGcorr_RSD.*STGcorr)],'VariableNames', ['RunName','Time','Sample','Elapse', STGTE_elements]);


end
