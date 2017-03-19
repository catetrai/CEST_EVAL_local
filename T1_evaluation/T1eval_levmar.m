function [T1info, T1map, popt] = T1eval_levmar(T1_stack,mapflag,ROInumber,P,Segment,StartValues)
% ** T1eval_levmar(T1_stack,mapflag,ROInumber,P,Segment,StartValues)
% mapflag: if 1 fits full given image (inside Segment)
% ROInumber: number of ROIs fitted
% P : P struct


if nargin<6
    StartValues = false;
    sprintf('No StartValues were given')
else
    sprintf('Using the following starting parameters:')
    StartValues
end

try
    TI=P.SEQ.w;
catch
    error('No inversion times given!');
end

if (nargin<5 && mapflag==1)
    sprintf('No Segment was given to function T1eval')
    Segment = make_Segment(T1_stack,'free');
    assignin('base','Segment',Segment);
end

if size(T1_stack,4) ~= numel(TI)
    error('Number of inversion times in P.SEQ.w is different from size of T1_stack')
end

% -------------------------------------------------------------------------
if mapflag
    
    [popt, P] = FIT_3D(T1_stack,P,Segment,StartValues);
    T1map=popt(:,:,:,1);
    
    figure,
    for ii=1:size(T1_stack,3)
        subplot(1,size(T1_stack,3),ii);
        imagesc(T1map(:,:,ii),[0 5000]);
        axis image
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        colormap jet;
        if (ii == size(T1_stack,3))
            colorbar
        end
        title('T1 map')
    end
    
else
    T1map=1;
    popt=0;
end;

% ROI eval
if (ROInumber~=0)
    % define ROI_def for all ROIs
%     if (size(T1_stack,3)>1)
%         prompt = {'Select slice for ROI analysis:'};
%         dlg_title = 'Input';
%         answer = inputdlg(prompt,dlg_title,1,{'1'});    
%         selected_slice=abs(str2num(answer{1}));
%     else
        selected_slice=1;
%     end
    
    ROIS = ROItool(squeeze(T1_stack(:,:,selected_slice,end)),ROInumber,'free');
    
     % preview plot
    if mapflag
        if (ndims(T1map)==2)
            im = T1map;
        else
            im = T1map(:,:,selected_slice);
        end
    else
        im=squeeze(T1_stack(:,:,selected_slice,end));
    end
    
    im=double(im);
    im(isnan(im))=0;
    
    figure;
    subplot(1,2,1);
    imagesc(im, [0.1*mean(im(im~=0)) 2*mean(im(im~=0))]), axis image;
    title('ROI definition')

    for ii=1:ROInumber
        % get ROI_def for current ROI
        ROI_def=ROIS{ii}.ROI_def;
        % fit every pixel in ROI
        [ROI_popt, P] = FIT_3D(T1_stack(:,:,selected_slice,:),P,ROI_def,StartValues);
        T1mapROI=squeeze(ROI_popt(:,:,:,1));
        
        % write values (mean,std,...) for current ROI
        T1info{ii}.mean = mean(T1mapROI(ROI_def==1));
        T1info{ii}.std = std(T1mapROI(ROI_def==1));
        T1info{ii}.img = T1mapROI;
        T1info{ii}.ROI_def = ROI_def;
        
        % calculate mean values of fit-result for preview plot
        for kk=1:size(ROI_popt,4);
            temp=ROI_popt(:,:,1,kk);
            ROI_popt_mean(kk)=mean(temp(ROI_def==1));
        end
        
        % calculate mean signal and std for all TI in current ROI
        for jj=1:numel(TI)
            buffer=squeeze(T1_stack(:,:,selected_slice,jj));
            ROI_mean(jj)= mean(buffer(ROI_def == 1));
            ROI_std(jj)= std(double(buffer(ROI_def == 1)));
        end;
       
        % calculate y_fit values for current ROI
        [f] = fitmodelfunc_ANA(ROI_popt_mean,P);
        
        % save all ROI results into S-struct
        S{ii}.y = ROI_mean;
        S{ii}.y_std = ROI_std;
        S{ii}.popt = ROI_popt_mean;
        S{ii}.T1 = T1info{ii}.mean;
        S{ii}.dT1 = T1info{ii}.std;
        S{ii}.y_fit = f(TI);
        
        TI_inter=[min(TI):1:max(TI)];
    
        subplot(1,2,1);
        hold on;
        contour(T1info{ii}.ROI_def,1,'m-','LineWidth',2);
        [yy xx]=find(T1info{ii}.ROI_def);
        t_h=text(fix(max(xx)),fix(min(yy)),'# ');
        set(t_h,'String',sprintf('ROI # %d',ii),'BackgroundColor',[1 1 1])
        leg1{ii}=sprintf('# ROI %d',ii);

        subplot(1,2,2)
        hold on
        leg2{ii}=sprintf('ROI%i: T1= %4.0f +- %4.0f ms',ii,S{ii}.T1,S{ii}.dT1);
        %leg{numel(leg)+1}=sprintf('ROI%i: data',ii);
        plot(TI_inter,f(TI_inter),'Color',cl(ii,ROInumber));
        hold on;
        h2(ii)=plot(TI,ROI_mean,':o','Color',cl(ii,ROInumber));
        
    end;
%    subplot(1,2,1)
%    legend(leg1,'location','SouthEast'); 
   subplot(1,2,2)
   legend(h2,leg2,'location','East'); 
   box on 
    
     
else
    T1info = 1;
end