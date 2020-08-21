clear all
%Script for opening and reading all raw agilent data
%Creates files including the collated raw counts data, standard deviations,
%peak integration times, number of cycles.
dirlist = dir('C:\Users\mdumo\OneDrive - University of St Andrews\Agilent\Matt\Full agilent data');
dirlist(ismember( {dirlist.name}, {'.', '..'})) = [];
Folders = dirlist([dirlist.isdir]);

counter=1;
for j=1:size(Folders,1)
    
    template=fullfile(Folders(j).folder,Folders(j).name);
   Batchlogloc=fullfile(template,'BatchLog.csv');  
if exist(Batchlogloc, 'file')~=2
    continue
end
 optsbatch = detectImportOptions(Batchlogloc); 
    optsbatch = setvaropts(optsbatch,'Acq_Date_Time','InputFormat','MM/dd/uuuu hh:mm:ss aa');
    batchlog=readtable(Batchlogloc,optsbatch);
    Samplelist=batchlog.SampleName(ismember(batchlog.AcquisitionResult,'Pass'));
    if isempty(Samplelist)
 continue
    end 
    runname=dirlist(j).name;   
  filelist=cell2mat(batchlog.FileName(ismember(batchlog.AcquisitionResult,'Pass') & ~ismember(batchlog.FileName,'-')));
  subFolders=string(filelist(:,end-8:end));
      
for i = 1:size(subFolders,1)

      file=dir(fullfile(template,subFolders(i),'*1.csv'));
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

%%
%Put all of the run tables together
  if counter == 1
      %Create table from the first sequence
     raw_CPS_T=run_CPS_t; 
     raw_SDs_T=run_SDs_t;
     raw_N_T=run_N_t;
     raw_intTime_T=run_intTime_t;   
  else
      %Further sequences added concatination, adding new elements if needed
     a=setdiff(run_CPS_t.Properties.VariableNames,raw_CPS_T.Properties.VariableNames);
     b=setdiff(raw_CPS_T.Properties.VariableNames,run_CPS_t.Properties.VariableNames);  

     t2_data = [run_CPS_t array2table(nan(height(run_CPS_t), numel(b)), 'VariableNames', b)];
     T2_data = [raw_CPS_T array2table(nan(height(raw_CPS_T), numel(a)), 'VariableNames', a)];
     
       t2_SDs=[run_SDs_t array2table(nan(height(run_SDs_t), numel(b)), 'VariableNames', b)];
       T2_SDs=[raw_SDs_T array2table(nan(height(raw_SDs_T), numel(a)), 'VariableNames', a)];
        
       t2_N=[run_N_t array2table(nan(height(run_N_t), numel(b)), 'VariableNames', b)];
       T2_N=[raw_N_T array2table(nan(height(raw_N_T), numel(a)), 'VariableNames', a)];
        
       t2_intTime=[run_intTime_t array2table(nan(height(run_intTime_t), numel(b)), 'VariableNames', b)];
       T2_intTime=[raw_intTime_T array2table(nan(height(raw_intTime_T), numel(a)), 'VariableNames', a)];
     
       raw_CPS_T=[T2_data;t2_data];
       raw_SDs_T=[T2_SDs;t2_SDs];
       raw_N_T=[T2_N;t2_N];
       raw_intTime_T=[T2_intTime;t2_intTime];
         
  end
    counter=counter+1;
end %%%%%%%%%%%%%% end of cycling through runs
%%

%load handel;
%player = audioplayer(y, Fs); 
%play(player);

save('Compiler_Ts.mat', 'raw_CPS_T', 'raw_N_T', 'raw_SDs_T', 'raw_intTime_T');
raw_CPS_T=sortrows(raw_CPS_T,{'Time'},{'ascend'});
raw_N_T=sortrows(raw_N_T,{'Time'},{'ascend'});
raw_intTime_T=sortrows(raw_intTime_T,{'Time'},{'ascend'});
raw_SDs_T=sortrows(raw_SDs_T,{'Time'},{'ascend'});
writetable(raw_CPS_T, 'raw_CPS_T.csv');
writetable(raw_SDs_T, 'raw_SDs_T.csv');
writetable(raw_N_T, 'raw_N_T.csv');
writetable(raw_intTime_T, 'raw_intTime_T.csv');
