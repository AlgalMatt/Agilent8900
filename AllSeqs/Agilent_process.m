% This script processes all data from the Agilent 8900. Agilent compiler.m
% is required.
% You will notice that the data is constantly converted from tables to
% matrices then back again to tables. This requires a lot more code than just
% working with matrices but it makes visual inspection of the data at each
% step much easier.

%03/08/20
%Make blank correction to nearest neighbour or linear not overall mean


clear all
close all
%Reads tables created by Agilent_conpiler.m
load Compiler_Ts.mat
Elements=raw_CPS_T.Properties.VariableNames(4:end); %Pulls full list of elements
runnames=unique(string(raw_CPS_T.RunName), 'stable'); %all sequence (run) names

%theoretical error vs real error
N_array=raw_N_T{:, Elements}; %number of cycles
CPS_array=raw_CPS_T{:, Elements}; %raw CPS data
realSE=raw_SDs_T{:, Elements}./raw_N_T{:, Elements}.^0.5;
realRSE=realSE./CPS_array;
realRSD=raw_SDs_T{:, Elements}./CPS_array; %relative standard deviation
realRSD_T=[raw_CPS_T(:,1:3), array2table(realRSD,'VariableNames', Elements)];

rawcounts = CPS_array.*raw_intTime_T{:, Elements}; %total counts
rawcounts(rawcounts<=0)=NaN;
theoSE=rawcounts.^0.5; %theoretical standard dev for single element analysis
theoRSE=theoSE./CPS_array; %theoretical relative standard dev
theoSE_T=[raw_CPS_T(:,1:3), array2table(theoSE,'VariableNames', Elements)];
theoRSE_T=[raw_CPS_T(:,1:3), array2table(theoRSE,'VariableNames', Elements)];

BLKROWS=contains(lower(raw_CPS_T.Sample), 'blk'); %blank locations



%Pre-allocate tables for computational efficiency
bmean_T=cell2table(cell(0,numel(raw_CPS_T.Properties.VariableNames)), ...
    'VariableNames', raw_CPS_T.Properties.VariableNames);
bmean_cali_T=bmean_T;
bcorr_T=cell2table(cell(0,numel(raw_CPS_T.Properties.VariableNames)+1), ...
    'VariableNames', [raw_CPS_T.Properties.VariableNames(1:3), 'Elapse', raw_CPS_T.Properties.VariableNames(4:end)]);
bcorr_SD_T=bcorr_T;
Ca43r_T=bcorr_T; Neff_T=bcorr_T;
Ca43r_RSD_T=bcorr_T; Ca43r_RSE_T=bcorr_T;
Ca43r_theoRSD_T=bcorr_T;
NeffRSD_T=bcorr_T; Neff_T=bcorr_T;

%Published STGTE values
STGTE_ratios=[20.711, 86.1, 8.09, 5.382, 5.382, 42.382, 58.634,...
    1.488, 0.14, 1.277, 5.629, 159.567];
STGTE_elements={'Li7', 'B11', 'Na23', 'Mg24', 'Mg25', 'Al27', 'Mn55', ...
    'Sr88', 'Cd111', 'Ba138', 'Nd146', 'U238'};

%pre-allocation of STGTE tables for computational efficiency
STGcorr_T=cell2table(cell(0,numel(STGTE_elements)+4), ...
    'VariableNames', [NeffRSD_T.Properties.VariableNames(1:4), STGTE_elements]);
STGcorr_RSD_T=STGcorr_T; STGcorr_theoRSD_T=STGcorr_T;
ec_raw_T=STGcorr_T; ec_blk_T=STGcorr_T; ec_Ca_T= STGcorr_T; ec_STG_T=STGcorr_T;

for i=1:numel(runnames) %start of cyling through each run
    runrows = find(ismember(raw_CPS_T.RunName, runnames(i)));
    %pull data from given run
    CPS_t=raw_CPS_T(runrows,:);
    elapse=CPS_t.Time-CPS_t.Time(1); %Time since beginning of run
    elapse(elapse>duration(60,00,00))=duration(00,00,00); %max limit on elapse. Stops duplicate analyses.
    info=[table2cell(CPS_t(:,1:3)), num2cell(elapse)];
    realSD_t=raw_SDs_T(runrows,Elements);
    run_rawcounts = rawcounts(runrows,:);
    run_N=N_array(runrows, :);
    
    %Find all the consistency standards, blanks and STGTE.
    blkrows=find(BLKROWS(runrows));
    STGTErows=find(contains(lower(CPS_t.Sample),'stgte'));
    
    
    
    %% blank corrections
    if isempty(blkrows)
        continue
    end
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
    bmean_SD=(nansum(blks_SD.^2,1).^0.5)./nansum(~isnan(blks_SD), 1);
    %place mean blank data back into a table
    bmean_t=cell2table([table2cell(CPS_t(blkrows(1),1:3)), num2cell(bmean)],...
        'VariableNames', ['RunName','Time','Sample', Elements]);
    
    
    %{
%Blank correction by mean blank
bcorr=CPS_t{:,Elements}-bmean;
%error propagation (could index ~BLKrows, but would need to re-jig sample names in table also)
bcorr_realSD=(bmean_SD.^2+realSD_t{:,Elements}.^2).^0.5;
bcorr_RSD=bcorr_realSD./bcorr;
%tabulate
bcorr_t=cell2table([info,num2cell(bcorr)],'VariableNames', ['RunName','Time','Sample','Elapse', Elements]);
bcorr_RSD_t=cell2table([info,num2cell(bcorr_RSD)],'VariableNames', ['RunName','Time','Sample','Elapse', Elements]);
        %}
        
        
        
        %pre-allocate tables for computational efficiency
        bcorr=zeros(size(CPS_t{:,Elements}));
        bcorr_SD=bcorr;
        
        secs=seconds(elapse);
        for j=1:numel(secs)
            [val, idx]=sort(blkrows_ol-j, 'ComparisonMethod', 'abs');
            %If there are not two STGTE bracketing the sample then just use the
            %closest single BLK as a point correction.
            
            bcorr(j,:)=CPS_t{j,Elements}-CPS_t{blkrows_ol(idx(1)),Elements};
            bcorr_SD(j,:)=(realSD_t{j,Elements}.^2+realSD_t{blkrows_ol(idx(1)),Elements}.^2).^0.5;
            %Otherwise use linear regression
            if ~all(val>=0) && ~all(val<=0)
                uidx=blkrows_ol(idx(find(val>0, 1)));
                lidx=blkrows_ol(idx(find(val<0, 1)));
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
        bcorr_RSD_t=cell2table([info,num2cell(bcorr_RSD)],'VariableNames', ['RunName','Time','Sample','Elapse', Elements]);
        
        %% Me/43Ca ratios for the standards
        Ca43_loc=find(ismember(Elements,'Ca43'));
        Ca43r=bcorr./bcorr(:,Ca43_loc);
        Ca43r_t=cell2table([info,num2cell(Ca43r)],'VariableNames', ['RunName','Time','Sample', 'Elapse', Elements]);
        Ca43r_RSD=(bcorr_RSD.^2+bcorr_RSD(:, Ca43_loc).^2).^0.5;
        Ca43r_RSD(isinf(Ca43r_RSD))=NaN;
        
        %% Me/Ca counting stats
        blkcorr_cpc = bcorr.*raw_intTime_T{runrows,Elements};
        %From John and Adkins (2010)
        Neff = (blkcorr_cpc.*blkcorr_cpc(:,Ca43_loc))./(blkcorr_cpc+blkcorr_cpc(:,Ca43_loc));
        %remove negatives and zeros
        Neff(Neff<=0)=NaN;
        %NeffRSD=(Neff.^0.5)./Neff;
        NeffRSD=1./(Neff.^0.5);
        Neff_t=cell2table([info,num2cell(Neff)],'VariableNames', ['RunName','Time','Sample', 'Elapse', Elements]);
        NeffRSD_t=cell2table([info,num2cell(NeffRSD)],'VariableNames', ['RunName','Time','Sample', 'Elapse', Elements]);
        %realSD_blankcorr=(bstd.^2+realSD_array.^2).^0.5;
        %Ca43Ratio_RSD=((realSD_blankcorr./blkcorr_array).^2+(realSD_blankcorr./blkcorr_array).^2).^0.5;
        Ca43r_RSE_t=cell2table([info,num2cell((Ca43r_RSD.*Ca43r./(run_N.^0.5))./Ca43r)],'VariableNames', ['RunName','Time','Sample', 'Elapse', Elements]);
        Ca43r_RSD_t=cell2table([info,num2cell(Ca43r_RSD)],'VariableNames', ['RunName','Time','Sample', 'Elapse', Elements]);
        %Ca43r_theoRSD_t=cell2table([info,num2cell(Ca43r_theoRSD)],'VariableNames', ['RunName','Time','Sample', Elements]);
        
        
        
        %% Put all run data into tables for easier indexing
        bmean_T=[bmean_T;bmean_t];
        bcorr_T=[bcorr_T;bcorr_t];
        bcorr_SD_T=[bcorr_SD_T; bcorr_SD_t];
        Ca43r_T=[Ca43r_T;Ca43r_t];
        Neff_T=[Neff_T;Neff_t];
        NeffRSD_T=[NeffRSD_T;NeffRSD_t];
        Ca43r_RSD_T=[Ca43r_RSD_T;Ca43r_RSD_t];
        Ca43r_RSE_T=[Ca43r_RSE_T;Ca43r_RSE_t];
        %Ca43r_theoRSD_T=[Ca43r_theoRSD_T;Ca43r_theoRSD_t];
        
        
        %% standard calibration
        if isempty(STGTErows)
            continue
        end
        %find locations of STGTE
        STGTErows=find(contains(lower(Ca43r_t.Sample),'stgte'));
        
        %Discount any 0.5STGTE
        if ~all(contains(lower(Ca43r_t.Sample(STGTErows)), '0.5stgte'))
            STGTErows(contains(lower(Ca43r_t.Sample(STGTErows)), '0.5stgte'))=[];
        end
        
        %remove anomalously low counts (blockages or exhaustion of sample)
        Ca43r_stg=Ca43r_t{:,STGTE_elements};
        Ca43r_RSDstg=Ca43r_RSD_t{:,STGTE_elements};
        Ca43r_stg(bcorr_t{:,STGTE_elements}<bmean_t{:,STGTE_elements}*3 | Ca43r_stg<0)=nan;
        Ca43r_RSDstg(isnan(Ca43r_stg))=nan;
        
        %pre-allocate tables for computational efficiency
        STGcorr=zeros(size(Ca43r_stg));
        STGcorr_RSD=STGcorr;
        
        secs=seconds(Ca43r_t{:, 'Elapse'});

        for j=1:numel(secs)
            [val, idx]=sort(STGTErows-j, 'ComparisonMethod', 'abs');
            %If there are not two STGTE bracketing the sample then just use the
            %closest single STGTE as a point correction.
            STGcorr(j,:)=Ca43r_stg(j,:)./Ca43r_stg(STGTErows(idx(1)),:).*STGTE_ratios;
            STGcorr_RSD(j,:)=(Ca43r_RSDstg(j,:).^2+Ca43r_RSDstg(STGTErows(idx(1)),:).^2).^0.5;                 
            %Otherwise use linear regression
            if ~all(val>=0) && ~all(val<=0)
                uidx=STGTErows(idx(find(val>0, 1)));
                lidx=STGTErows(idx(find(val<=0, 1)));
                c=Ca43r_stg(lidx, :);
                m=(Ca43r_stg(uidx, :)-c)./(secs(uidx)-secs(lidx));
                STGcorr(j,:)=Ca43r_stg(j,:)./(m.*(secs(j)-secs(lidx))+c).*STGTE_ratios;
                
                %error propagation (regress)
                c=Ca43r_RSDstg(lidx, :);
                m=(Ca43r_RSDstg(uidx,:)-c)./(secs(uidx)-secs(lidx));
                
                STGcorr_RSD(j,:)=(Ca43r_RSDstg(j,:).^2+(m.*(secs(j)-secs(lidx))+c).^2).^0.5;
            end
        end
        
        %% contribution of errors
        runrealRSDSTG=realRSD_T{runrows,STGTE_elements};
        nb_rows=~contains(lower(CPS_t.Sample), 'blk'); %discount blanks
        info_nb=info(nb_rows, :);
        %relative error contributions at each stage proportional to total error
        ec_raw=runrealRSDSTG(nb_rows,:)./STGcorr_RSD(nb_rows,:);
        ec_blk=(bcorr_RSD_t{nb_rows,STGTE_elements}-runrealRSDSTG(nb_rows,:))./STGcorr_RSD(nb_rows,:);
        ec_Ca=(Ca43r_RSD_t{nb_rows,STGTE_elements}-bcorr_RSD_t{nb_rows,STGTE_elements})./STGcorr_RSD(nb_rows,:);
        ec_STG=(STGcorr_RSD(nb_rows,:)-Ca43r_RSD_t{nb_rows,STGTE_elements})./STGcorr_RSD(nb_rows,:);
        
        ec_raw_t=cell2table([info_nb,num2cell(ec_raw)],'VariableNames', ['RunName','Time','Sample', 'Elapse', STGTE_elements]);
        ec_blk_t=cell2table([info_nb,num2cell(ec_blk)],'VariableNames', ['RunName','Time','Sample', 'Elapse', STGTE_elements]);
        ec_Ca_t=cell2table([info_nb,num2cell(ec_Ca)],'VariableNames', ['RunName','Time','Sample', 'Elapse', STGTE_elements]);
        ec_STG_t=cell2table([info_nb,num2cell(ec_STG)],'VariableNames', ['RunName','Time','Sample', 'Elapse', STGTE_elements]);
        
        %Remove outliers
        STGcorr_RSD(STGcorr_RSD>1E3)=nan;
        
         STGcorr_t=cell2table([info,num2cell(STGcorr)],'VariableNames', ['RunName','Time','Sample', 'Elapse', STGTE_elements]);
        STGcorr_RSD_t=cell2table([info,num2cell(STGcorr_RSD)],'VariableNames', ['RunName','Time','Sample','Elapse', STGTE_elements]);
        STGcorr_SD_t=cell2table([info,num2cell(STGcorr_RSD.*STGcorr)],'VariableNames', ['RunName','Time','Sample','Elapse', STGTE_elements]);
        
        %remove blanks
        %STGcorr_t= STGcorr_t(nb_rows, :);
        %STGcorr_RSD_t=STGcorr_RSD_t(nb_rows, :);
        %STGcorr_SD_t=STGcorr_SD_t(nb_rows, :);
        
        %Join the run table with the rest of the runs.
        STGcorr_T=[STGcorr_T;STGcorr_t];
        STGcorr_RSD_T=[STGcorr_RSD_T;STGcorr_RSD_t];
        STGcorr_SD_T=[STGcorr_RSD_T;STGcorr_SD_t];
        ec_raw_T=[ec_raw_T; ec_raw_t];
        ec_blk_T=[ec_blk_T; ec_blk_t];
        ec_Ca_T= [ec_Ca_T; ec_Ca_t];
        ec_STG_T= [ec_STG_T; ec_STG_t];
        
        
        %% Calibrate blank to STGTE
        %used to observe normalised blank variability through time
        
        bmean_cali=bmean./nanmean(CPS_t{STGTErows,Elements},1) ;
        info_b=info(contains(lower(CPS_t.Sample), 'blk'), :);
        bmean_cali_t=cell2table([info_b(1,1:3),num2cell(bmean_cali)],'VariableNames', ['RunName','Time','Sample', Elements]);
        %add to a table
        bmean_cali_T=[bmean_cali_T;bmean_cali_t];
        
end


  bol=bmean_cali_T{:, Elements}>nanmean(bmean_cali_T{:, Elements},1)...
      +2.5*nanstd(bmean_cali_T{:, Elements},1) | bmean_cali_T{:, ...
      Elements}<nanmean(bmean_cali_T{:, Elements},1)-...
      2.5*nanstd(bmean_cali_T{:, Elements},1);
    while ~isempty(find(bol, 1))
        bmean_cali_T{:, Elements}(bol)=NaN;
         bol=bmean_cali_T{:, Elements}>nanmean(bmean_cali_T{:, Elements},1)...
      +2.5*nanstd(bmean_cali_T{:, Elements},1) | bmean_cali_T{:, ...
      Elements}<nanmean(bmean_cali_T{:, Elements},1)-...
      2.5*nanstd(bmean_cali_T{:, Elements},1);
    end







%Extract individual consistency standards
CS1_T=STGcorr_T(contains(lower(STGcorr_T.Sample), 'cs1'),:);
CS2_T=STGcorr_T(contains(lower(STGcorr_T.Sample), 'cs2'),:);
CS3_T=STGcorr_T(contains(lower(STGcorr_T.Sample), 'cs3'),:);
std8301f_T=STGcorr_T(contains(lower(STGcorr_T.Sample), '8301f'),:);
std8301c_T=STGcorr_T(contains(lower(STGcorr_T.Sample), '8301c'),:);

%Remove outliers
for i=1:numel(STGTE_elements)
       
    CS1_array(:,i)=CS1_T{:,STGTE_elements(i)};
    outs=CS1_array(:,i)>nanmean(CS1_array(:,i))+2.5*nanstd(CS1_array(:,i)) | CS1_array(:,i)<nanmean(CS1_array(:,i))-2.5*nanstd(CS1_array(:,i))| isinf(CS1_array(:,i));
    while ~all(~outs)
    CS1_array(outs, i)=NaN;
    outs=CS1_array(:,i)>nanmean(CS1_array(:,i))+2.5*nanstd(CS1_array(:,i)) | CS1_array(:,i)<nanmean(CS1_array(:,i))-2.5*nanstd(CS1_array(:,i))| isinf(CS1_array(:,i));
    end
       
      CS2_array(:,i)=CS2_T{:,STGTE_elements(i)};
    outs=CS2_array(:,i)>nanmean(CS2_array(:,i))+2.5*nanstd(CS2_array(:,i)) | CS2_array(:,i)<nanmean(CS2_array(:,i))-2.5*nanstd(CS2_array(:,i))| isinf(CS2_array(:,i));
    while ~all(~outs)
    CS2_array(outs, i)=NaN;
    outs=CS2_array(:,i)>nanmean(CS2_array(:,i))+2.5*nanstd(CS2_array(:,i)) | CS2_array(:,i)<nanmean(CS2_array(:,i))-2.5*nanstd(CS2_array(:,i))| isinf(CS2_array(:,i));
    end
    
      CS3_array(:,i)=CS3_T{:,STGTE_elements(i)};
    outs=CS3_array(:,i)>nanmean(CS3_array(:,i))+2.5*nanstd(CS3_array(:,i)) | CS3_array(:,i)<nanmean(CS3_array(:,i))-2.5*nanstd(CS3_array(:,i))| isinf(CS3_array(:,i));
    while ~all(~outs)
    CS3_array(outs, i)=NaN;
    outs=CS3_array(:,i)>nanmean(CS3_array(:,i))+2.5*nanstd(CS3_array(:,i)) | CS3_array(:,i)<nanmean(CS3_array(:,i))-2.5*nanstd(CS3_array(:,i))| isinf(CS3_array(:,i));
    end
    
      std8301f_array(:,i)=std8301f_T{:,STGTE_elements(i)};
    outs=std8301f_array(:,i)>nanmean(std8301f_array(:,i))+2.5*nanstd(std8301f_array(:,i)) | std8301f_array(:,i)<nanmean(std8301f_array(:,i))-2.5*nanstd(std8301f_array(:,i))| isinf(std8301f_array(:,i));
    while ~all(~outs)
    std8301f_array(outs, i)=NaN;
    outs=std8301f_array(:,i)>nanmean(std8301f_array(:,i))+2.5*nanstd(std8301f_array(:,i)) | std8301f_array(:,i)<nanmean(std8301f_array(:,i))-2.5*nanstd(std8301f_array(:,i))| isinf(std8301f_array(:,i));
    end
  
      std8301c_array(:,i)=std8301c_T{:,STGTE_elements(i)};
    outs=std8301c_array(:,i)>nanmean(std8301c_array(:,i))+2.5*nanstd(std8301c_array(:,i)) | std8301c_array(:,i)<nanmean(std8301c_array(:,i))-2.5*nanstd(std8301c_array(:,i))| isinf(std8301c_array(:,i));
    while ~all(~outs)
    std8301c_array(outs, i)=NaN;
    outs=std8301c_array(:,i)>nanmean(std8301c_array(:,i))+2.5*nanstd(std8301c_array(:,i)) | std8301c_array(:,i)<nanmean(std8301c_array(:,i))-2.5*nanstd(std8301c_array(:,i))| isinf(std8301c_array(:,i));
    end
    

    %Number of each consistency standards
    CS1_n(i)=length(CS1_array(~isnan(CS1_array(:,i)), i));
    CS2_n(i)=length(CS2_array(~isnan(CS2_array(:,i)), i));
    CS3_n(i)=length(CS3_array(~isnan(CS3_array(:,i)), i));
    std8301f_n(i)=length(std8301f_array(~isnan(std8301f_array(:,i)), i));
    std8301c_n(i)=length(std8301c_array(~isnan(std8301c_array(:,i)), i));
    
    %Mean & stddev of each consistency standards
    CS1_mean(i)=nanmean(CS1_array(:,i));
    CS1_SD(i)=nanstd(CS1_array(:,i));
    
    CS2_mean(i)=nanmean(CS2_array(:,i));
    CS2_SD(i)=nanstd(CS2_array(:,i));
    
    CS3_mean(i)=nanmean(CS3_array(:,i));
    CS3_SD(i)=nanstd(CS3_array(:,i));
    
    std8301f_mean(i)=nanmean(std8301f_array(:,i));
    std8301f_SD(i)=nanstd(std8301f_array(:,i));

    std8301c_mean(i)=nanmean(std8301c_array(:,i));
    std8301c_SD(i)=nanstd(std8301c_array(:,i));
    
end

%Standard data without outliers
CS1_T=[CS1_T(:,1:4) cell2table(num2cell(CS1_array),'VariableNames', STGTE_elements)];
CS2_T=[CS2_T(:,1:4) cell2table(num2cell(CS2_array),'VariableNames', STGTE_elements)];
CS3_T=[CS3_T(:,1:4) cell2table(num2cell(CS3_array),'VariableNames', STGTE_elements)];
std8301f_T=[std8301f_T(:,1:4) cell2table(num2cell(std8301f_array),'VariableNames', STGTE_elements)];
std8301c_T=[std8301c_T(:,1:4) cell2table(num2cell(std8301c_array),'VariableNames', STGTE_elements)];

opts = detectImportOptions('NISTRM8301.xlsx');
opts.RowNamesRange='A2:A12';
opts.Sheet='8301cdat';
i8301c=readtable('NISTRM8301.xlsx', opts);
opts.Sheet='8301fdat';
i8301f=readtable('NISTRM8301.xlsx', opts);
opts.Sheet='8301crsd';
i8301c_rsd=readtable('NISTRM8301.xlsx', opts);
opts.Sheet='8301frsd';
i8301f_rsd=readtable('NISTRM8301.xlsx', opts);

save('Agilentprocess.mat')
save('Consistency.mat', 'CS1_T', 'CS2_T', 'CS3_T', 'std8301f_T', 'std8301c_T', 'bmean_cali_T');
save('Interlab.mat', 'i8301c', 'i8301c_rsd', 'i8301f', 'i8301f_rsd');









