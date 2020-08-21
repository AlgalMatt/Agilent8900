close all
clear all
load Agilentprocess.mat
load agilentsim.mat
load Interlab.mat
%theoretical error vs real error
N_array=raw_N_T{:, Elements}; %number of cycles
CPS_array=raw_CPS_T{:, Elements}; %raw CPS data
realSE=raw_SDs_T{:, Elements}./(raw_N_T{:, Elements}.^0.5);
realRSE=realSE./CPS_array;
realRSD=raw_SDs_T{:, Elements}./CPS_array; %relative standard deviation
realRSD_T=[raw_CPS_T(:,1:3), array2table(realRSD,'VariableNames', Elements)];

CPS_array(CPS_array<=0)=NaN; % lambda
CPC = CPS_array.*raw_intTime_T{:, Elements}; %counts per cycle
totcounts=CPC.*N_array; %total counts

theoSD=CPS_array.^0.5; %Standard dev for Poisson distribution lambda^0.5
theoSE=(CPS_array./raw_intTime_T{:, Elements}).^0.5; %SE for Poisson distribution (lambda/n)^0.5
theoRSD=(CPC.^0.5)./CPC;
theoRSE=theoSE./CPS_array; %theoretical relative standard err

%Extract individual consistency standards
CS1_indx=contains(lower(raw_CPS_T.Sample), 'cs1');
CS2_indx=contains(lower(raw_CPS_T.Sample), 'cs2');
CS3_indx=contains(lower(raw_CPS_T.Sample), 'cs3');



CPC_CS1=CPC(CS1_indx, :);
CPC_CS2=CPC(CS2_indx, :);
CPC_CS3=CPC(CS3_indx, :);

for i=1:numel(Elements)
    
    outs=CPC_CS1(:, i)>nanmean(CPC_CS1(:, i))+2.5*nanstd(CPC_CS1(:, i)) | CPC_CS1(:, i)<nanmean(CPC_CS1(:, i))-2.5*nanstd(CPC_CS1(:, i))| isinf(CPC_CS1(:, i));
    while ~all(~outs)
        CPC_CS1(outs, i)=NaN;
        outs=CPC_CS1(:, i)>nanmean(CPC_CS1(:, i))+2.5*nanstd(CPC_CS1(:, i)) | CPC_CS1(:, i)<nanmean(CPC_CS1(:, i))-2.5*nanstd(CPC_CS1(:, i))| isinf(CPC_CS1(:, i));
    end
    
    outs=CPC_CS2(:, i)>nanmean(CPC_CS2(:, i))+2.5*nanstd(CPC_CS2(:, i)) | CPC_CS2(:, i)<nanmean(CPC_CS2(:, i))-2.5*nanstd(CPC_CS2(:, i))| isinf(CPC_CS2(:, i));
    while ~all(~outs)
        CPC_CS2(outs, i)=NaN;
        outs=CPC_CS2(:, i)>nanmean(CPC_CS2(:, i))+2.5*nanstd(CPC_CS2(:, i)) | CPC_CS2(:, i)<nanmean(CPC_CS2(:, i))-2.5*nanstd(CPC_CS2(:, i))| isinf(CPC_CS2(:, i));
    end
    
    
    outs=CPC_CS3(:, i)>nanmean(CPC_CS3(:, i))+2.5*nanstd(CPC_CS3(:, i)) | CPC_CS3(:, i)<nanmean(CPC_CS3(:, i))-2.5*nanstd(CPC_CS3(:, i))| isinf(CPC_CS3(:, i));
    while ~all(~outs)
        CPC_CS3(outs, i)=NaN;
        outs=CPC_CS3(:, i)>nanmean(CPC_CS3(:, i))+2.5*nanstd(CPC_CS3(:, i)) | CPC_CS3(:, i)<nanmean(CPC_CS3(:, i))-2.5*nanstd(CPC_CS3(:, i))| isinf(CPC_CS3(:, i));
    end
end

%mean CS counts
CS_CPCmean=nanmean(CPC_CS1);
CS_CPCmean(2,:)=nanmean(CPC_CS2);
CS_CPCmean(3,:)=nanmean(CPC_CS3);


optimalint=15000./CS_CPCmean.*raw_intTime_T{999:1001, 4:end};



BLKROWS_r=contains(lower(Ca43r_T.Sample), 'blk');

%% error error plot (elements)
for i=1:numel(Elements)
    fig0=figure;
    plot(theoRSD(~BLKROWS, i), realRSD(~BLKROWS, i), 'ob');
    hold on
    plot(theoRSD(BLKROWS, i), realRSD(BLKROWS, i), '*r');
    plot([0 2], [0 2], '-k', 'LineWidth',2)
    xlabel('theoretical error')
    ylabel('real error')
    title(Elements{i})
    ax=gca;
    ax.YLim=[0 1];
    ax.XLim=ax.YLim;
    ax.YLim=[0 max(realRSD(:, i))*1.1];
    ax.XLim=[0 max(theoRSD(:, i))*1.1];
    shadedErrorBar(theotot,mu,sdmu*2,'lineProps','g'); %shaded errors
    legend('samples', 'blanks', 'Location', 'best')
    hold off
    print(fig0,strcat('error_errorfigs/Me_',Elements{i}, '_SD'),'-dpng')
    close all
end



%An attempt to show error bars that are representative of the Poisson
%distribution, not normal distribution
for i=1:numel(Elements)
    fig011=figure;
    plot(theoRSD(~BLKROWS, i), realRSD(~BLKROWS, i), 'ob');
    hold on
    plot(theoRSD(BLKROWS, i), realRSD(BLKROWS, i), '*r');
    plot([0 2], [0 2], '-k', 'LineWidth',2)
    xlabel('theoretical error')
    ylabel('real error')
    title(Elements{i})
    ax=gca;
    ax.YLim=[0 1];
    ax.XLim=ax.YLim;
    ax.YLim=[0 max(realRSD(:, i))*1.1];
    ax.XLim=[0 max(theoRSD(:, i))*1.1];
    shadedErrorBar(theotot,mu,sdmu_ee,'lineProps','g'); %shaded errors
    legend('samples', 'blanks', 'Location', 'best')
    hold off
    print(fig011,strcat('error_errorfigs/Me_',Elements{i}, '_SDP'),'-dpng')
    close all
end







%% error error plot (ratios)
for i=1:numel(Elements)
    elementplot=Elements{i};
    Fig1pt5=figure;
    x=NeffRSD_T{~BLKROWS_r, elementplot};
    y=Ca43r_RSD_T{~BLKROWS_r, elementplot};
    plot(x, y, 'ob');
    hold on
    plot([0 1], [0 1], '-k', 'LineWidth',2)
    xlabel('theoretical error (1/Neff^{0.5})')
    ylabel('real error (propagated)')
    title([Elements{i},'/Ca43'])
    ax=gca;
    if ~isnan(max(y))
    ax.YLim=[0 max(y(y<nanmean(y(y<nanmean(y)+nanstd(y)))+nanstd(y(y<nanmean(y)+nanstd(y)))))];
    ax.XLim=[0 max([max(x(y<nanmean(y(y<nanmean(y)+nanstd(y)))+nanstd(y(y<nanmean(y)+nanstd(y))))), 0.2])];
    end
    hold off
    print(Fig1pt5,strcat('error_errorfigs/CaMe_',Elements{i}),'-dpng')
    close all
end





%% total counts per run vs error
for i=1:numel(Elements)
    elementplot=Elements{i};
    %%%%Raw counts plot
    slopex=min(log10(CPC(:,i))):0.1:max(log10(CPC(:,i)));
    slopey=-0.5*slopex;
    Fig2=figure('units','normalized','outerposition',[0 0 1 1]);
    a=  plot(log10(CPC(~BLKROWS,i)), log10(realRSD(~BLKROWS,i)), 'ob');
    hold on
    b= plot(log10(CPC(BLKROWS,i)), log10(realRSD(BLKROWS,i)), '*r');
    c=  plot(slopex, slopey, 'k', 'LineWidth', 2) ;
    shadedErrorBar(log10(simtotcounts/5),log10(mu),[log10(mu+2*sdmu)-log10(mu);log10(mu)-log10(mu-sdmu)],'lineProps','g');
    x=log10(CPC(:,i)); y=log10(realRSD(:,i));
    idx=isnan(x) | isnan(y);
    xy=sortrows([x(~idx), y(~idx)]);
    pval=[];
    for j=1:numel(x(~idx))
        mdl=fitlm(xy(j:end, 1), xy(j:end, 2), 'linear');
        pval(j)=mdl.Coefficients.pValue(2);
    end
    inflect=xy(find(pval>0.05, 1), 1);
    d= plot(inflect, interp1(slopex, slopey, inflect), '*c', 'MarkerSize', 20, 'LineWidth', 1.5);
    p1=polyfit(x(~idx), y(~idx), 2);
    polyy=polyval(p1, slopex);
    e=plot(slopex, polyy, 'm--', 'LineWidth', 2);
    
    
    idx=isnan(x) | isnan(y) |  x>3.5;
    p2=polyfit(x(~idx), y(~idx), 1);
    polyy=polyval(p2, slopex);
    f=plot(slopex, polyy, '--', 'color',[0.5 0.5 0.5], 'LineWidth', 2);
    
    xlabel('log(tot counts per cycle)')
    ylabel('log (RSE)')
    title(elementplot)
    if ~isempty(nanmean(y(x>6))) && ~isnan(nanmean(y(x>6)))
        g=yline(nanmean(y(x>6)), ':k', 'LineWidth', 2);
        legend([a b c d e f g], 'samples', 'blanks', 'y=0.5x','inflection pt (pvalue > 0.05)',...
            '2 deg poly fit', 'linear fit (x < 3.5)' , 'mean (x > 6)');
        
    else
        legend([a b c e f], 'samples', 'blanks', 'inflection pt (pvalue > 0.05)',...
            '2 deg poly fit', 'linear fit (x < 3.5)');
    end
    
    xlim([min(slopex), max(slopex)])
    ylim([min(log10(realRSD(:,i))), max(log10(realRSD(:,i)))])
    hold off
    print(Fig2,['Count_stats\Elements\',elementplot, '_CountsPlot'],'-dpng')
    close all
end




%on average there is a 14555 variance


normrnd(100,1.4455, 1, 1000)


%cause of Ca scatter at high counts
hiidx=realRSD(:, 7)>0.01 & log10(CPC(:,7))>5;
loidx=realRSD(:, 7)<0.01 & log10(CPC(:,7))>5;
lorsdCa=raw_CPS_T(loidx, :);
hirsdCa=raw_CPS_T(hiidx, :);

idx=log10(CPC(:,7))>5;
Cadat=realRSD(idx, 7);
xdat=CPS_array(idx, :);

plot(Cadat, datenum(raw_CPS_T.Time(idx)), 'o')
xlim([0 0.1])
plot(Cadat, hour(raw_CPS_T.Time(idx)), 'o')
xlim([0 0.1])
seqorder=[];
for i=1:numel(runnames)
    seqorder=[seqorder; [1:numel(raw_CPS_T.Sample(raw_CPS_T.RunName==runnames(i)))]'];
end

plot(Cadat, seqorder(idx), 'o')
xlim([0 0.1])

runnames2=unique(raw_CPS_T{idx, 'RunName'});
for j=28
    i=7;
    idx2=ismember(raw_CPS_T{:, 'RunName'}, runnames2(j));
    elementplot=Elements{i};
    %%%%Raw counts plot
    slopex=min(log10(CPC(:,i))):0.1:max(log10(CPC(:,i)));
    slopey=-0.5*slopex;
    Fig2=figure('units','normalized','outerposition',[0 0 1 1]);
    a=  plot(log10(CPC(~BLKROWS,i)), log10(realRSD(~BLKROWS,i)), 'ob');
    hold on
    b= plot(log10(CPC(BLKROWS,i)), log10(realRSD(BLKROWS,i)), '*r');
    c=plot(slopex, slopey, 'k', 'LineWidth', 2) ;
    d=yline(nanmean(y(x>6)), ':k', 'LineWidth', 2);
    e=plot(log10(CPC(idx2,i)), log10(realRSD(idx2,i)), '*k', 'LineWidth', 3);
    xlabel('log(tot counts)')
    ylabel('log (RSE)')
    title(runnames2(j))
    hold off
end

for i=1:numel(Elements)
    figure;
    plot(Cadat, xdat(:, i), 'o')
    hold on
    title(Elements{i})
    xlim([0 0.06])
    ylim([0 nanmean(xdat(:, i))+2*nanstd(xdat(:, i))])
    hold off
end


%% Neff vs error
for i=1:numel(Elements)
    elementplot=Elements{i};
    
    slopex=min(log10(Neff_T{~BLKROWS_r,elementplot})):max(log10(Neff_T{~BLKROWS_r,elementplot}));
    slopey=-0.5*slopex;
    Fig1=figure('units','normalized','outerposition',[0 0 1 1]);
    plot(log10(Neff_T{~BLKROWS_r,elementplot}), log10(Ca43r_RSD_T{~BLKROWS_r,elementplot}), 'o')
    hold on
    plot(slopex, slopey, 'r', 'LineWidth', 3)
    xlabel('log(Neff)')
    ylabel('log (RSE)')
    title([elementplot, '/Ca43'])
    hold off
    print(Fig1,['Count_stats\Ratios\',elementplot, '_NeffPlot'],'-dpng')
    close all
end




%% time vs normalised blank counts as a fraction of STGTE (1mM)

for i=1:numel(Elements)
    elementplot=Elements{i};
    fig8=figure;
    plot(bmean_cali_T.Time,bmean_cali_T{:,i+3},'o');
    title(Elements{i})
    ax=gca;
    if diff(ax.YLim)>1
        set(ax, 'YScale', 'log')
    end
    xlabel('Time')
    ylabel('Blank CPS/STGTE CPS')
    print(fig8,['timefigs\normblanks\',elementplot, '_blanks'],'-dpng')
    close all
end



%% Standards over time




for i=1:numel(STGTE_elements)
    %CS1 time vs Me/Ca
    fig3=figure;
    plot(CS1_T.Time, CS1_array(:,i), 'o')
    xlabel('Time')
    ylabel(strcat(STGTE_elements(i),'/Ca43 (\mumol/mol)'))
    yline(CS1_mean(i), 'r');
    yline(CS1_mean(i)+CS1_SD(i)*2, '-.r');
    yline(CS1_mean(i)-CS1_SD(i)*2, '-.r');
    %ylim([min([CS1_mean(i)-nanstd(CS1_array(:,i))*3 CS1_mean(i)-3]) max([CS1_mean(i)+nanstd(CS1_array(:,i))*3 CS1_mean(i)+3])]);
    ylim([CS1_mean(i)-CS1_SD(i)*3 CS1_mean(i)+CS1_SD(i)*3])
    title(['mean=' num2str(CS1_mean(i), '%.2f') ', 2RSD(%)=' num2str(CS1_SD(i)*2/CS1_mean(i)*100, '%.2f') ', n=' num2str(CS1_n(i))]);
    
    %CS2 time vs Me/Ca
    fig4=figure;
    plot(CS2_T.Time, CS2_array(:,i), 'o')
    xlabel('Time')
    ylabel(strcat(STGTE_elements(i),'/Ca43 (\mumol/mol)'))
    yline(CS2_mean(i), 'r');
    yline(CS2_mean(i)+CS2_SD(i)*2, '-.r');
    yline(CS2_mean(i)-CS2_SD(i)*2, '-.r');
    %ylim([min([CS2_mean(i)-nanstd(CS2_array(:,i))*3 CS2_mean(i)-3]) max([CS2_mean(i)+nanstd(CS2_array(:,i))*3 CS2_mean(i)+3])]);
    ylim([CS2_mean(i)-CS2_SD(i)*3 CS2_mean(i)+CS2_SD(i)*3])
    title(['mean=' num2str(CS2_mean(i), '%.2f') ', 2RSD(%)=' num2str(CS2_SD(i)*2/CS2_mean(i)*100, '%.2f') ', n=' num2str(CS2_n(i))]);
    
    %CS3 time vs Me/Ca
    fig5=figure;
    plot(CS3_T.Time, CS3_array(:,i), 'o')
    xlabel('Time')
    ylabel(strcat(STGTE_elements(i),'/Ca43 (\mumol/mol)'))
    yline(CS3_mean(i), 'r');
    yline(CS3_mean(i)+CS3_SD(i)*2, '-.r');
    yline(CS3_mean(i)-CS3_SD(i)*2, '-.r');
    %ylim([min([CS3_mean(i)-nanstd(CS3_array(:,i))*3 CS3_mean(i)-3]) max([CS3_mean(i)+nanstd(CS3_array(:,i))*3 CS3_mean(i)+3])]);
    ylim([CS3_mean(i)-CS3_SD(i)*3 CS3_mean(i)+CS3_SD(i)*3])
    title(['mean=' num2str(CS3_mean(i), '%.2f') ', 2RSD(%)=' num2str(CS3_SD(i)*2/CS3_mean(i)*100, '%.2f') ', n=' num2str(CS3_n(i))]);
    
    print(fig3,strcat('timefigs/CS1/CS1_',STGTE_elements{i}),'-dpng')
    print(fig4,strcat('timefigs/CS2/CS2_',STGTE_elements{i}),'-dpng')
    print(fig5,strcat('timefigs/CS3/CS3_',STGTE_elements{i}),'-dpng')
end
for i=1:numel(STGTE_elements)
    %8301F time vs Me/Ca
    fig6=figure('units','normalized','outerposition',[0 0 1 1]);
    el=STGTE_elements{i};
    plot(datenum(std8301f_T.Time), std8301f_T{:,STGTE_elements(i)}, 'ob')    
   hold on
    datetick('x',1)
    xlabel('Time')
    ylabel([el,'/Ca43 (\mumol/mol)'])
    yline(std8301f_mean(i), 'b');
    yline(std8301f_mean(i)+std8301f_SD(i), '-.b');
    yline(std8301f_mean(i)-std8301f_SD(i), '-.b');
    %ylim([min([CS3_mean(i)-nanstd(CS3_array(:,i))*3 CS3_mean(i)-3]) max([CS3_mean(i)+nanstd(CS3_array(:,i))*3 CS3_mean(i)+3])]);
    title(['mean=' num2str(std8301f_mean(i), '%.2f') ', 2RSD(%)=' num2str(std8301f_SD(i)*2/std8301f_mean(i)*100, ...
        '%.2f') ', n=' num2str(std8301f_n(i))]);
    %interlab stuff
    if STGTE_elements{i}=="Mg25"
       el='Mg24'; 
    end
    yline(i8301f{'Mean', el}, 'r');
     ax=gca;
     rectangle('Position', [ax.XLim(1), i8301f{'Mean', el}-i8301f{'sd', el}, ...
         ax.XLim(2)-ax.XLim(1), 2*i8301f{'sd', el}], 'FaceColor',[0.5 0 0, 0.2], ...
         'EdgeColor', [0, 0, 0, 0]);
    errorbar(repmat(ax.XLim(1)+10, 7,1), i8301f{1:7, el}, i8301f_rsd{1:7, el}.*i8301f{1:7, el},...
        'sr', 'MarkerFaceColor', 'r')
    text(repmat(ax.XLim(1)+10, 7,1), i8301f{1:7, el},...
        [" NIST"; " UoB"; " GEOMAR"; " LSCE"; " Ox"; " Soton"; " StA"], 'Color', 'r')
    hold off
     print(fig6,strcat('timefigs/std8301f/8301f_',STGTE_elements{i}),'-dpng')
    close all
end



%{
for i=1:numel(STGTE_elements)
    %8301c time vs Me/Ca
    fig6=figure('units','normalized','outerposition',[0 0 1 1]);
    el=STGTE_elements{i};
    plot(datenum(std8301c_T.Time), std8301c_T{:,STGTE_elements(i)}, 'ob')    
   hold on
    datetick('x',1)
    xlabel('Time')
    ylabel([el,'/Ca43 (\mumol/mol)'])
    yline(std8301c_mean(i), 'b');
    yline(std8301c_mean(i)+std8301c_SD(i), '-.b');
    yline(std8301c_mean(i)-std8301c_SD(i), '-.b');
    %ylim([min([CS3_mean(i)-nanstd(CS3_array(:,i))*3 CS3_mean(i)-3]) max([CS3_mean(i)+nanstd(CS3_array(:,i))*3 CS3_mean(i)+3])]);
    title(['mean=' num2str(std8301c_mean(i), '%.2f') ', 2RSD(%)=' num2str(std8301c_SD(i)*2/std8301c_mean(i)*100, ...
        '%.2f') ', n=' num2str(std8301c_n(i))]);
    %interlab stuff
    if STGTE_elements{i}=="Mg25"
       el='Mg24'; 
    end
    yline(i8301c{'Mean', el}, 'r');
     ax=gca;
     rectangle('Position', [ax.XLim(1), i8301c{'Mean', el}-i8301c{'sd', el}, ...
         ax.XLim(2)-ax.XLim(1), 2*i8301c{'sd', el}], 'FaceColor',[0.5 0 0, 0.2], ...
         'EdgeColor', [0, 0, 0, 0]);
    errorbar(repmat(ax.XLim(1)+10, 7,1), i8301c{1:7, el}, i8301c_rsd{1:7, el}.*i8301c{1:7, el},...
        'sr', 'MarkerFaceColor', 'r')
    text(repmat(ax.XLim(1)+10, 7,1), i8301c{1:7, el},...
        [' Lab1'; ' Lab2'; ' Lab3'; ' Lab4'; ' Lab5'; ' Lab6'; ' Lab7'], 'Color', 'r')
    hold off
     print(fig6,strcat('timefigs/std8301c/8301c_',STGTE_elements{i}),'-dpng')
    close all
end
%}
%% error contributions

for i=1:numel(STGTE_elements)
    %Box plots of error contributions
    fig7=figure;
    boxplot([ec_raw_T{:, i+3}, ec_blk_T{:, i+3}, ec_Ca_T{:, i+3}, ec_STG_T{:, i+3}], ...
        'Labels',{'sample analyte counts','blank', 'sample Ca counts', 'STGTE correction'})
    ylabel('Fractional error contribution')
    title(STGTE_elements{i})
    print(fig7,strcat('error_contributions/',STGTE_elements{i}),'-dpng')
    
    close all
end




%% Run drift

stg_rows=find(contains(lower(Ca43r_T.Sample),'stgte'));

%Discount any 0.5STGTE
if ~all(contains(lower(Ca43r_T.Sample(stg_rows)), '0.5stgte'))
    stg_rows(contains(lower(Ca43r_T.Sample(stg_rows)), '0.5stgte'))=[];
end
dur=[duration(0, 5, 0), duration(0, 10, 0), duration(0, 15, 0),...
    duration(0, 20, 0), duration(0, 25, 0), duration(0, 30, 0),...
    duration(0, 40, 0), duration(0, 50, 0), duration(1, 00, 0)];

  


for i=1:numel(STGTE_elements)
        stgT=Ca43r_T(stg_rows, :);
    yrsd=Ca43r_RSD_T(stg_rows, :);
    %remove outliers
    outs=stgT{:, STGTE_elements(i)}>nanmean(stgT{:, STGTE_elements(i)})+3*nanstd(stgT{:, STGTE_elements(i)}) | stgT{:, STGTE_elements(i)}<nanmean(stgT{:, STGTE_elements(i)})-3*nanstd(stgT{:, STGTE_elements(i)})| isinf(stgT{:, STGTE_elements(i)});
    while ~all(~outs)
        stgT(outs, :)=[];
        yrsd(outs, :)=[];
        outs=stgT{:, STGTE_elements(i)}>nanmean(stgT{:, STGTE_elements(i)})+3*nanstd(stgT{:, STGTE_elements(i)}) | stgT{:, STGTE_elements(i)}<nanmean(stgT{:, STGTE_elements(i)})-3*nanstd(stgT{:, STGTE_elements(i)})| isinf(stgT{:, STGTE_elements(i)}) | stgT{:, STGTE_elements(i)}<0;
    end
    
    stgT=sortrows(stgT, 'RunName'); yrsd=sortrows(yrsd, 'RunName');
    y2sd=(yrsd{:, STGTE_elements(i)}.*stgT{:, STGTE_elements(i)})*2;
    runs=unique(stgT.RunName);
    
    yall=[];
    %These plots show the change in STGTE measurements of element/Ca ratios
    %over each run relative to the first STGTE measurement made during that
    %run. This is accompanied by an annotation of the mean error (2sd) of
    %individual STGTE measurements of that ratio (i.e. the error within any
    %given single measurement of STGTE). Measurements of STGTE that fall outside
    %of this threshold should be considered as significantly drifted and
    %therefore warrant regular inclusions of STGTE throughout a run to correct
    %for this drift.
    
    decnum=0; %number of runs that show an overall decrease (not necessarily significant)
    driftnum=0; %number of runs that show a significant drift
    fig8=figure('units','normalized','outerposition',[0 0 1 1]);
    for j=1:numel(runs)
        yrel=stgT{ismember(stgT.RunName, runs(j)), STGTE_elements(i)}- stgT{find(ismember(stgT.RunName, runs(j)), 1), STGTE_elements(i)};
        x=stgT{ismember(stgT.RunName, runs(j)), 'Elapse'};
        yall=[yall;yrel];
        if yrel(end)<0
            decnum=decnum+1;
        end
        p1=plot(x, yrel, 'o-');
        if any(yrel>nanmean(y2sd)) || any(yrel<-nanmean(y2sd))
           driftnum=driftnum+1;
            t=text(x(end), yrel(end), ['  ', datestr(stgT{find(ismember(stgT.RunName, runs(j)), 1), 'Time'}, 'dd/mm/yy')]);
            t.Color=p1.Color;   
        end
        hold on
    end
    ax=gca;
    lims(i, :)=ax.YLim;
    yline(nanmean(y2sd), '--k', ['Average 2\sigma for STGTE ', STGTE_elements{i}, '/Ca43']);
    yline(-nanmean(y2sd), '--k', ['Average 2\sigma for STGTE ', STGTE_elements{i}, '/Ca43']);
    title([STGTE_elements{i}]);
   annotation('textbox',[0.68 0.6 0.3 0.3],'String',[num2str(decnum/numel(runs)*100, 2) ,...
        '% of runs showed decrease. ', num2str(driftnum/numel(runs)*100, 2), '% of runs drifted'], 'FitBoxToText','on');
    xlabel('Time since start of run')
    ylabel(['\Delta',STGTE_elements{i},'/Ca43'])
    hold off
     print(fig8,['Drift/',STGTE_elements{i}],'-dpng')
    %errorbar(datenum(stgT{:, 4}), yall,  y2sd,  'o');
    close Figure 1
    
end 


%The above only tells us that having a single STGTE at the beginning of the
%run will result in a bias. 
%How often are STGTE needed to properly correct?
%Would linear regression rather than point correction work better?




for i=1:numel(STGTE_elements)  
    stgT=Ca43r_T(stg_rows, :);
    yrsd=Ca43r_RSD_T(stg_rows, :);
    %remove outliers
    outs=stgT{:, STGTE_elements(i)}>nanmean(stgT{:, STGTE_elements(i)})+3*nanstd(stgT{:, STGTE_elements(i)}) | stgT{:, STGTE_elements(i)}<nanmean(stgT{:, STGTE_elements(i)})-3*nanstd(stgT{:, STGTE_elements(i)})| isinf(stgT{:, STGTE_elements(i)});
    while ~all(~outs)
        stgT(outs, :)=[];
        yrsd(outs, :)=[];
        outs=stgT{:, STGTE_elements(i)}>nanmean(stgT{:, STGTE_elements(i)})+3*nanstd(stgT{:, STGTE_elements(i)}) | stgT{:, STGTE_elements(i)}<nanmean(stgT{:, STGTE_elements(i)})-3*nanstd(stgT{:, STGTE_elements(i)})| isinf(stgT{:, STGTE_elements(i)}) | stgT{:, STGTE_elements(i)}<0;
    end
    
    stgT=sortrows(stgT, 'RunName'); yrsd=sortrows(yrsd, 'RunName');
    y2sd=(yrsd{:, STGTE_elements(i)}.*stgT{:, STGTE_elements(i)})*2;
    runs=unique(stgT.RunName);
    
    
    decnum=0;
    driftnum=0;
    fig9=figure('units','normalized','outerposition',[0 0 1 1]);
    for j=1:numel(runs)
        y=stgT{ismember(stgT.RunName, runs(j)), STGTE_elements(i)};
        x=stgT{ismember(stgT.RunName, runs(j)), 'Elapse'};
        
        xsec=seconds(x);
        stgint=[];
        stggap=80*10; %interval betweem STGTE corrections through run (seconds)
        for k=0:floor(xsec(end)/stggap)
            [val, idx]=sort(xsec-stggap*k, 'ComparisonMethod', 'abs');
            
            stgint(k+1)=idx(1);
            
        end
        stgint(end+1)=numel(x);
        yresid=[];
        for k=1:numel(x)
            [val, idx]=sort(stgint-k, 'ComparisonMethod', 'abs');
            if val(1)~=0
                m=(y(stgint(idx(find(val>0, 1))))-y(stgint(idx(find(val<0, 1)))))/(xsec(stgint(idx(find(val>0, 1))))-xsec(stgint(idx(find(val<0, 1)))));
                c=y(stgint(idx(find(val<0, 1))));
                yresid(k)=y(k)-(m*(xsec(k)-xsec(stgint(idx(find(val<0, 1)))))+c);
            else
                yresid(k)=0;
            end
            
        end
        yresid2(1:numel(yresid), j)=yresid;
        p1=plot(x, yresid, 'o-');
        %if any(yrel>nanmean(y2sd)) || any(yrel<-nanmean(y2sd))
        %   driftnum=driftnum+1;
        %  t=text(x(end), yrel(end), ['  ', datestr(stgT{find(ismember(stgT.RunName, runs(j)), 1), 'Time'}, 'dd/mm/yy')]);
        % t.Color=p1.Color;
        %end
        hold on
    end
    yline(nanmean(y2sd), '--k', ['Average 2\sigma for STGTE ', STGTE_elements{i}, '/Ca43']);
    yline(-nanmean(y2sd), '--k', ['Average 2\sigma for STGTE ', STGTE_elements{i}, '/Ca43']);
    yline(nanmean(y2sd)/2, ':k', ['Average 1\sigma for STGTE ', STGTE_elements{i}, '/Ca43']);
    yline(-nanmean(y2sd)/2, ':k', ['Average 1\sigma for STGTE ', STGTE_elements{i}, '/Ca43']);
    title([STGTE_elements{i}]);
    % annotation('textbox',[0.68 0.6 0.3 0.3],'String',[num2str(decnum/numel(runs)*100, 2) ,...
    %    '% of runs showed decrease. ', num2str(driftnum/numel(runs)*100, 2), '% of runs drifted'], 'FitBoxToText','on');
    xlabel('Time since start of run')
    ylabel(['\Delta',STGTE_elements{i},'/Ca43'])
    ylim(lims(i,:))
    hold off
    print(fig9,['Drift/Residuals_',STGTE_elements{i}],'-dpng')
    %errorbar(datenum(stgT{:, 4}), yall,  y2sd,  'o');
    close Figure 1 
end



%% Magnesium



MgCS1=bcorr_T{contains(lower(bcorr_T.Sample), 'cs1'), 'Mg24'}./...
    bcorr_T{contains(lower(bcorr_T.Sample), 'cs1'), 'Mg25'};
Mg8301f=bcorr_T{contains(lower(bcorr_T.Sample), '8301f'), 'Mg24'}./...
    bcorr_T{contains(lower(bcorr_T.Sample), '8301f'), 'Mg25'};
MgSTG=bcorr_T{contains(lower(bcorr_T.Sample), '1stgte'), 'Mg24'}./...
    bcorr_T{contains(lower(bcorr_T.Sample), '1stgte'), 'Mg25'};
MgCS3=bcorr_T{contains(lower(bcorr_T.Sample), 'cs3'), 'Mg24'}./...
    bcorr_T{contains(lower(bcorr_T.Sample), 'cs3'), 'Mg25'};

 plot(bcorr_T{contains(lower(bcorr_T.Sample), '8301f'), 'Time'}, Mg8301f, 'or')
hold on
plot(bcorr_T{contains(lower(bcorr_T.Sample), '1stgte'), 'Time'}, MgSTG, 'ob')
plot(bcorr_T{contains(lower(bcorr_T.Sample), 'cs1'), 'Time'}, MgCS1, 'ok')
plot(bcorr_T{contains(lower(bcorr_T.Sample), 'cs3'), 'Time'}, MgCS3, 'og')
   xlabel('Time')
    ylabel(['Mg24/Mg25'])
legend('CS1 (Mg/Ca = 1.31)', '8301f (Mg/Ca = 2.59)', 'STGTE (Mg/Ca = 5.38)', 'CS3 (Mg/Ca = 7.58)')
ylim([6 10])
yline(7.9);
hold off













