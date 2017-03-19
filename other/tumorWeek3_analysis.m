%% PET NORMALIZATION


%% Read PMOD '.pixelDump' files
direc = dir('/Users/caterina/Desktop/CEST/WISC_backup/Caterina/PMOD/rats_tumor/PIXELDUMPS');
direc(cell2mat({direc.isdir})) = [];

for ixfile = 1:length(direc)
    ixrat = direc(ixfile).name(10);
    ixslice = find( find(~cellfun('isempty', regexp(({direc.name})', ['rat',ixrat]))) == ixfile );
    fid = fopen(['/Users/caterina/Desktop/CEST/WISC_backup/Caterina/PMOD/rats_tumor/PIXELDUMPS/',direc(ixfile).name]);
    d = textscan(fid, '%s%s%s%f%d%d%d%f%s%s%s', 'Headerlines', 3);
    istumor = strcmp(d{3},'3');
    ispet = cellfun('isempty', regexp(d{1},'MTRrex'));
    week3.(['rat',ixrat]).pet.t{ixslice} = d{8}(istumor & ispet);
    week3.(['rat',ixrat]).pet.h{ixslice} = d{8}(~istumor & ispet);
    week3.(['rat',ixrat]).cest.t{ixslice} = d{8}(istumor & ~ispet);
    week3.(['rat',ixrat]).cest.h{ixslice} = d{8}(~istumor & ~ispet);
    fclose(fid);
end
%% Example rat: look at histograms of CEST and PET for tumor vs. healthy
nrat = 'rat6';
figure;
for ixs = 1:size(week3.(nrat).pet.t,2)
    subplot(2,4,ixs);
    histogram(week3.(nrat).pet.t{ixs}, 10, 'Normalization', 'probability');
    hold on; histogram(week3.(nrat).pet.h{ixs}, 10, 'Normalization', 'probability');
    legend('tumor','healthy')
    title([nrat, ' PET'])
    subplot(2,4,ixs+4);
    histogram(week3.(nrat).cest.t{ixs}, 10, 'Normalization', 'probability');
    hold on; histogram(week3.(nrat).cest.h{ixs}, 10, 'Normalization', 'probability');
    legend('tumor','healthy')
    title([nrat, ' CEST'])
end

%% Boxplots of example rat, to look at variability between slices (CEST and PET)
nrat = 'rat6';
figure;
for ixs = 1:size(week3.(nrat).pet.t,2)
    
    subplot(2,4,ixs);
    g_pet = [repmat('t',[size(week3.(nrat).pet.t{ixs}),1]); repmat('h',[size(week3.(nrat).pet.h{ixs}),1])];
    boxplot([week3.(nrat).pet.t{ixs}; week3.(nrat).pet.h{ixs}], g_pet, 'labels', {'tumor', 'control'});
    alphabonf = 0.05/length(direc);
    [H,P] = kstest2(week3.(nrat).pet.t{ixs}, week3.(nrat).pet.h{ixs}, 'alpha', alphabonf);
    if H && P<0.001
%         title(sprintf('%s PET (p<0.001)', nrat))
        title(sprintf('Slice %d (p<0.001)', ixs))
    else
%         title([nrat, ' PET (n.s.)'])
        title(sprintf('Slice %d (n.s.)', ixs))
    end
    set(gca,'fontsize',16)
    
    week3.(nrat).cest.t{ixs} = week3.(nrat).cest.t{ixs}/6.5534e+04;
    week3.(nrat).cest.h{ixs} = week3.(nrat).cest.h{ixs}/6.5534e+04;
    
    subplot(2,4,ixs+4);
    g_cest = [repmat('t',[size(week3.(nrat).cest.t{ixs}),1]); repmat('h',[size(week3.(nrat).cest.h{ixs}),1])];
    boxplot([week3.(nrat).cest.t{ixs}; week3.(nrat).cest.h{ixs}], g_cest, 'labels', {'tumor', 'control'});
    [H,P] = kstest2(week3.(nrat).cest.t{ixs}, week3.(nrat).cest.h{ixs}, 'alpha', alphabonf);
    if H && P<001
        title(sprintf('Slice %d (p<0.001)', ixs))
    else
        title(sprintf('Slice %d (n.s.)', ixs))
    end
    set(gca,'fontsize',16)
end


%% ttest
i=0; t=[]; h=[];
for n=[1:4, 6]
    i=i+1;
    ratname=sprintf('rat%d',n);
    t(i)=mean(stats_cestunits.(ratname).cest.t(:,1));
    h(i)=mean(stats_cestunits.(ratname).cest.h(:,1));
end
    

%%
figure;errorbar(stats.(ratname).cest.t(:,1), stats.(ratname).cest.t(:,2), 'o')
hold on; errorbar(stats.(ratname).cest.h(:,1), stats.(ratname).cest.h(:,2), 'o')
legend('tumor','control')
title('amide CEST (MTR_{Rex} at 3.5 ppm)')
xlabel('slice #')


%% FIGURES
ratname='rat4';

% CEST
figure;errorbar(stats.(ratname).cest.t(:,1), stats.(ratname).cest.t(:,2), 'o')
hold on; errorbar(stats.(ratname).cest.h(:,1), stats.(ratname).cest.h(:,2), 'o')
legend('tumor','control')
title('amide CEST (MTR_{Rex} at 3.5 ppm)')
xlabel('slice #')

% PET
figure;errorbar(stats_bwnorm.(ratname).pet.t(:,1), stats_bwnorm.(ratname).pet.t(:,2), 'o')
hold on; errorbar(stats_bwnorm.(ratname).pet.h(:,1), stats_bwnorm.(ratname).pet.h(:,2), 'o')
legend('tumor','control')
title('[^{18}F]FMISO PET (kBq/cc)')
xlabel('slice #')

%% CEST-PET dependency
figure;
for n=[1:4, 6]
    ratname=sprintf('rat%d',n);
    plot(mean(stats.(ratname).cest.t(:,1)), mean(stats.(ratname).pet.t(:,1)), 'ro');
    hold on;
    plot(mean(stats.(ratname).cest.h(:,1)), mean(stats.(ratname).pet.h(:,1)), 'bo');
end
xlabel('amide CEST (MTR_{Rex})')
ylabel('{18}F]FMISO PET (%ID/cc)')

%% tumor/healthy ratio, SLICE-WISE
figure;
colors='rgbcky';
hrat=[];
for n=[1:4, 6]
    ratname=sprintf('rat%d',n);
    hrat(n)=plot(stats_cestunits.(ratname).cest.t(:,1)./stats_cestunits.(ratname).cest.h(:,1), ...
        stats_cestunits.(ratname).pet.t(:,1)./stats_cestunits.(ratname).pet.h(:,1),...
        'o', 'markerfacecolor', colors(n), 'MarkerEdgeColor', 'k', 'markersize', 8);
    hold on
end
xlabel('CEST (MTR_{Rex})')
ylabel('PET (SUV)')
% xlim([0 2]); ylim([0 2]);
title('Tumor/control ratio')
plot([0 2], [1 1], '--k')
plot([1 1], [0 2], '--k')
legend(hrat([1:4,6]), 'rat1', 'rat2', 'rat3', 'rat4', 'rat5')
set(gca, 'fontsize', 20)
axis equal
axis square

%% tumor/healthy ratio, SINGLE-ANIMAL (w/propagated errorbars)
ratioError = @(A,B,sd_A,sd_B) sqrt((1/B * sd_A)^2 + (-A*B^(-2) * sd_B)^2);

for n=1:5
    if n==5
        ratname=sprintf('rat%d',n+1);
    else
        ratname=sprintf('rat%d',n);
    end
    a_cest_week3(n) = mean(stats_cestunits.(ratname).cest.t(:,1));
    b_cest_week3(n) = mean(stats_cestunits.(ratname).cest.h(:,1));
    sd_a_cest_week3(n) = std(stats_cestunits.(ratname).cest.t(:,1));
    sd_b_cest_week3(n) = std(stats_cestunits.(ratname).cest.h(:,1));
    a_pet(n) = mean(stats_cestunits.(ratname).pet.t(:,1));
    b_pet(n) = mean(stats_cestunits.(ratname).pet.h(:,1));
    sd_a_pet(n) = std(stats_cestunits.(ratname).pet.t(:,1));
    sd_b_pet(n) = std(stats_cestunits.(ratname).pet.h(:,1));
    
    error_cest(n) = ratioError(a_cest_week3(n), b_cest_week3(n), sd_a_cest_week3(n), sd_b_cest_week3(n));
    error_pet(n) = ratioError(a_pet(n), b_pet(n), sd_a_pet(n), sd_b_pet(n));
end

%%
figure;
colors=cl(1:8,8); colors=colors(2:6,:);
hrat=[];
for i=1:5
    line([a_pet(i)/b_pet(i) - error_pet(i), a_pet(i)/b_pet(i) + error_pet(i)], [a_cest_week3(i)/b_cest_week3(i), a_cest_week3(i)/b_cest_week3(i)], 'color', 'k', 'linewidth', 1);
    hold on
    line([a_pet(i)/b_pet(i), a_pet(i)/b_pet(i)], [a_cest_week3(i)/b_cest_week3(i) - error_cest(i), a_cest_week3(i)/b_cest_week3(i) + error_cest(i)], 'color', 'k', 'linewidth', 1);
    hrat(i)=plot(a_pet(i)/b_pet(i), a_cest_week3(i)/b_cest_week3(i), 'o', 'markersize', 10, 'markerfacecolor', colors(i,:), 'MarkerEdgeColor', 'k');
end

plot([0 2], [1 1], '--k')
plot([1 1], [0 2], '--k')
legend(hrat, 'rat1', 'rat2', 'rat3', 'rat4', 'rat5', 'location','northwest')
ylabel('CEST (MTR_{Rex})')
xlabel('PET (SUV)')
title('Tumor/control ratio')
xlim([0 2]); ylim([0 2])
axis square
set(gca,'fontsize',20)


%% bar graph CEST
figure;
ratnames={'rat1', 'rat2', 'rat3', 'rat4', 'rat6'};
xaxis=[1 4 7 10 13];
for n=1:5
%     ratname=sprintf('rat%d',n);
    ht=bar(xaxis(n), mean(stats_cestunits.(ratnames{n}).cest.t(:,1)), 0.6, 'r', 'linewidth', 2);
    hold on
    errorbar(xaxis(n), mean(stats_cestunits.(ratnames{n}).cest.t(:,1)), 0, std(stats_cestunits.(ratnames{n}).cest.t(:,1)), 'k', 'linewidth', 2)
end
% errorbar(xaxis, mean(stats.(ratnames{n}).cest.t(:,1)), 0, std(stats.(ratnames{n}).cest.t(:,1)), 'k', 'linewidth', 2);

for n=1:5
%     ratname=sprintf('rat%d',n);
    hh=bar(xaxis(n)+1, mean(stats_cestunits.(ratnames{n}).cest.h(:,1)), 0.6, 'b', 'linewidth', 2);
    errorbar(xaxis(n)+1, mean(stats_cestunits.(ratnames{n}).cest.h(:,1)), 0, std(stats_cestunits.(ratnames{n}).cest.h(:,1)),'k', 'linewidth', 2)
end
title('amide CEST')
xlabel('animal #')
ylabel('MTR_{Rex}(3.5 ppm)')
set(gca,'xtick', xaxis+0.5, 'xticklabel', 1:5)
set(gca, 'fontsize', 20)
legend([ht hh], {'tumor', 'control'}, 'orientation','horizontal','location','northwest');

%% bar graph PET
figure;
ratname={'rat1', 'rat2', 'rat3', 'rat4', 'rat6'};
xaxis=[1 4 7 10 13];
for n=1:5
%     ratname=sprintf('rat%d',n);
    ht=bar(xaxis(n), mean(stats_cestunits.(ratname{n}).pet.t(:,1)), 0.6, 'r', 'linewidth', 2);
    hold on
    errorbar(xaxis(n), mean(stats_cestunits.(ratname{n}).pet.t(:,1)), 0, std(stats_cestunits.(ratname{n}).pet.t(:,1)), 'k', 'linewidth', 2)
end
for n=1:5
%     ratname=sprintf('rat%d',n);
    hh=bar(xaxis(n)+1, mean(stats_cestunits.(ratname{n}).pet.h(:,1)), 0.6, 'b', 'linewidth', 2);
    errorbar(xaxis(n)+1, mean(stats_cestunits.(ratname{n}).pet.h(:,1)), 0, std(stats_cestunits.(ratname{n}).pet.h(:,1)),'k', 'linewidth', 2)
end
title('[^{18}F]FMISO PET')
xlabel('animal #')
ylabel('SUV')
set(gca,'xtick', xaxis+0.5, 'xticklabel', 1:5)
set(gca, 'fontsize', 20)
legend([ht hh], {'tumor', 'control'}, 'orientation','horizontal','location','northwest');


%% TIME COURSE

figure;
xname={'week 1', 'week 2', 'week3'};
xaxis=[1 4 7];
for n=1:3
%     ratname=sprintf('rat%d',n);
    ht=bar(xaxis(n), a_cest(n), 0.6, 'r', 'linewidth', 2);
    hold on
    errorbar(xaxis(n), a_cest(n), 0, sd_a_cest(n), 'k', 'linewidth', 2)
end
for n=1:3
%     ratname=sprintf('rat%d',n);
    hh=bar(xaxis(n)+1, b_cest(n), 0.6, 'b', 'linewidth', 2);
    errorbar(xaxis(n)+1, b_cest(n), 0, sd_b_cest(n),'k', 'linewidth', 2)
end
title('amide CEST')
ylabel('MTR_{Rex}(3.5 ppm)')
set(gca,'xtick', xaxis+0.5, 'xticklabel', xname)
set(gca, 'fontsize', 20)
legend([ht hh], {'tumor', 'control'}, 'orientation','horizontal','location','northeast');


%% TIME COURSE RATIO
ratioError = @(A,B,sd_A,sd_B) sqrt((1/B * sd_A)^2 + (-A*B^(-2) * sd_B)^2);
for n=1:3
    error_cest(n) = ratioError(a_cest(n), b_cest(n), sd_a_cest(n), sd_b_cest(n));
end
figure;
errorbar(a_cest./b_cest, error_cest, 'ko', 'markersize',12, 'linewidth', 2.5)
set(gca,'xtick', 1:3, 'xticklabel', xname)
title('Tumor/control ratio')
set(gca, 'fontsize', 20)
hold on;
plot(get(gca,'xlim'), [1 1], '--r', 'linewidth', 1.8)