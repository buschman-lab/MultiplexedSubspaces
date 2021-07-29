function  [stPred,stTrue]=Deconvolve_GLM(dff,st,method,train_idx)
%Camden MacDowell - timeless
%set train index to 1 to train on everything
%also included nonglm methods in here since they are easy to run and the
%glm doesn't take long to fit. 

[n,z] = size(dff);

if nargin <4
    train_idx = false(n,1);
    train_idx(1:floor(n*0.7),:) = true;
end

if isempty(train_idx) %no comparisons 
    %fit kernel across all recordings
    win = 30; %size of the kernel window (def = 1 sec)
    dff = cat(1,dff(:));
    st = cat(1,st(:));    
    predictors = createRollingWindow(dff, win); %t-n:t-1
    response =  st(ceil(win/2):end-floor(win/2)); % get the middle timepoint in window  
    %fit kernel
    kernel = glmfit(predictors,response); %mean centering not needed here since you don't care about intercept
    kernel = kernel(2:end); %remove intercept  
    switch method
        case 'lucyrichardson' %req positive
            kernel = kernel-min(kernel);
            stPred = {deconvlucy(dff'-min(dff),kernel)};
        case 'lucyrichardson_gcamp' %easy to put in here instead of making another function
            stPred = {lucric(dff'-min(dff),0.95,1,30)};            
        case 'direct'
            stPred = {convn(padarray(dff',[0,floor(length(kernel)/2)],'replicate','both')',kernel,'valid')'};   
        otherwise
            error('unknown deconvolv method')
    end    
    stTrue = st';
    
else %run comparisons
    stPred = cell(1,z);
    for i = 1:z %loop through probes/recs
        fprintf('\n\t running kernel fits for probe/rec %d of %d', i,z);
        %fit kernel centered on predicted spiking        
        win = 31; %should be odd 
        if mod(win,2)==0; win = win+1;  end 
        dff_temp = dff(train_idx,i);
        st_temp = st(train_idx,i);
        predictors = createRollingWindow(dff_temp, win); %t-n:t-1
        response =  st_temp(ceil(win/2):end-floor(win/2)); % get the middle timepoint in window  
        %fit kernel
        kernel = glmfit(predictors,response); %mean centering not needed here since you don't care about intercept
        kernel = kernel(2:end); %remove intercept     
    
        stPred{i} = NaN(sum(train_idx==0),1);
        %test on all (including self)
        for j = 1:z %interal loop through probes/recs
            dff_test = dff(~train_idx,j)';
            switch method
                case 'lucyrichardson' %req positive
                    kernel = kernel-min(kernel);
                    stPred{i}(:,j) = deconvlucy(dff_test'-min(dff_test),kernel); 
                case 'lucyrichardson_gcamp' %use gcamp kernel
                    stPred{i}(:,j) = lucric(dff_test'-min(dff_test),0.95,1,30);
                case 'direct'
                    stPred{i}(:,j) = convn(padarray(dff_test,[0,floor(length(kernel)/2)],'replicate','both')',kernel,'valid');                
                otherwise
                    error('unknown deconvolv method')
            end      
            
        end %j loop        
        stTrue = st(~train_idx,:);
    end %i loop
    
    
end

end %function end




















        %
%         
%    temp = deconvlucy(dff_probe(:,cur_probe)-min(dff_probe(:,cur_probe)),kernel);
% %     temp = lucric(dff_probe(:,cur_probe)-min(dff_probe(:,cur_probe)),gp.d_gamma,gp.d_smooth,gp.d_kernel);
% %    temp = deconvreg(dff_probe(:,cur_probe),st_kernel);
% %    temp = convn(padarray(dff_probe(:,cur_probe)',[0,floor(length(st_kernel)/2)],'replicate','both')',st_kernel,'valid');           
%    
%         
%         
% %           st_kernel = kernel;   
%    st_kernel = st_kernel/max(st_kernel);
%    st_kernel = st_kernel-min(st_kernel);
%    temp = deconvlucy(dff_probe(:,cur_probe)-min(dff_probe(:,cur_probe)),st_kernel);
%         
%         
        
%         
%         
%         [~,~,netc,xic,aic] = train_narx_nn(dff(train_idx,i)',st(train_idx,i)');
%         
%         stPred = NaN(sum(train_idx==0),1);
%         %test on all (including self)
%         for j = 1:z %interal loop through probes/recs
%             dff_test = num2cell(dff(~train_idx,j)');
%             stPred(:,j) = cell2mat(netc(dff_test,xic,aic));                        
%         end %j loop
%         
%         stTrue = st(~train_idx,:)';
%     end %i loop
%     
%     %generate statistics
% %     stats = deconvolutionfitStats(stPred,stTrue); 
%     stats =[];
%     
%     %flip and return
%     stTrue = flipud(stTrue);
%     stPred = flipud(stPred);
%     
% end %comparison if/else
% 
% end %fuction 
%    
%    %predict imaging timing at timepoint t using preceding timepoints
%    win = round(opts.fps/2);
%    if mod(win,2)==0; win = win+1;  end 
%    predictors = createRollingWindow(dff_probe(1:end-1,cur_probe), win); %t-n:t-1
%    response = dff_probe(win+1:end,cur_probe);
%    auto_kernel = glmfit(predictors,response); %mean centering not needed here since you don't care about intercept
%    auto_kernel = auto_kernel(2:end); %remove intercept
%    %max normalize
%    auto_kernel = auto_kernel/max(auto_kernel);
%    
%    %regress out autocorr
%    dff_regressed = convn(padarray(dff_probe(:,cur_probe)',[0,floor(length(auto_kernel)/2)],'replicate','both')',auto_kernel,'valid');           
%    figure; plot(dff_regressed);  hold on; plot(dff_probe(:,cur_probe));
%    
%    %predict imaging timing at timepoint t using preceding ephys
%    predictors = createRollingWindow(st(1:end-1), win); %t-n:t-1
%    response = dff_regressed(win+1:end);
%    st_kernel = glmfit(predictors,response); %mean centering not needed here since you don't care about intercept
%    st_kernel = st_kernel(2:end); %remove intercept   
%    st_kernel = cat(1,-1*flipud(st_kernel),zeros(win-1,1));
%    st_kernel = st_kernel/max(st_kernel);
%     
%    %deconvolve
%    st_regressed = convn(padarray(st',[0,floor(length(st_kernel)/2)],'replicate','both')',flipud(st_kernel),'valid');           
%    
%    figure; plot(st_regressed); yyaxis right;  hold on; plot(dff_regressed);
%    
%    figure; plot(xcorr(dff_probe(:,cur_probe),st_regressed,200,'normalized'));
%    hold on; plot(xcorr(dff_regressed,st_regressed,200,'normalized'));
%    figure; histogram(st_regressed)
%    
%    %try the other direction
% %    temp = convn(padarray(dff_regressed',[0,floor(length(st_kernel)/2)],'replicate','both')',-1*(st_kernel),'valid');           
% %    figure; plot(xcorr(dff_probe(:,cur_probe),temp,200,'normalized'));
% %    hold on; plot(xcorr(dff_regressed,temp,200,'normalized'));
% %    hold on; plot(xcorr(dff_probe(:,cur_probe),st,200,'normalized'));
%    temp = deconv(dff_regressed,st_kernel);
%    figure; plot(temp); hold on; yyaxis right;  hold on; plot(st);
%    
%    parameter_class = 'general_params_corticaldynamics';
%    gp = loadobj(feval(parameter_class));   
%    lr = lucric(dff_probe(:,cur_probe),gp.d_gamma,gp.d_smooth,gp.d_kernel);
%    
%    figure; plot(temp); yyaxis right; hold on; plot(dff_regressed);
%       
%    figure; plot(lr); yyaxis right;  hold on; plot(st);
% %    figure; plot(xcorr(lr-nanmean(lr),200,'normalized'));
%    hold on; plot(xcorr(lr-nanmean(lr),st-nanmean(st),200,'normalized'));
%    hold on; plot(xcorr(temp-nanmean(temp),st-nanmean(st),200,'normalized'));
%    hold on; plot(xcorr(dff_probe(:,cur_probe)-nanmean(dff_probe(:,cur_probe)),st-nanmean(st),200,'normalized'));
%    
%    %close all
%    figure; hold on; 
%    hold on; plot(xcorr(dff_probe(:,cur_probe)-nanmean(dff_probe(:,cur_probe)),temp-nanmean(temp),200,'normalized'));
%    hold on; plot(xcorr(dff_probe(:,cur_probe)-nanmean(dff_probe(:,cur_probe)),st-nanmean(st),200,'normalized'));
%    hold on; plot(xcorr(temp-nanmean(temp),st-nanmean(st),200,'normalized'));
%    hold on; plot(xcorr(st-nanmean(st),st-nanmean(st),200,'normalized'));
% %    xcorr(dff_probe(:,cur_probe),st_regressed);
% %    dff_probe(:,cur_probe);
% %%  
% 
%    %predict st timing at timepoint t using following timepoints
%    win = round(opts.fps);
%    if mod(win,2)==0; win = win+1;  end 
%    predictors = createRollingWindow(dff_probe(2:end,cur_probe), win); %t-n:t-1
%    response = dff_probe(1:end-win,cur_probe);
%    auto_kernel = glmfit(predictors,response); %mean centering not needed here since you don't care about intercept
%    auto_kernel = auto_kernel(2:end); %remove intercept
%    %max normalize
% %    auto_kernel = auto_kernel/max(auto_kernel);
%    
%    %regress out autocorr
%    dff_regressed = convn(padarray(dff_probe(:,cur_probe)',[0,floor(length(auto_kernel)/2)],'replicate','both')',auto_kernel,'valid');           
%    figure; plot(dff_regressed);  hold on; plot(dff_probe(:,cur_probe));
%    
%    %predict imaging timing at timepoint t using preceding ephys
%    predictors = createRollingWindow(dff_probe(2:end,cur_probe), win); %t-n:t-1
%    response = st(1:end-win);
%    st_kernel = glmfit(predictors,response); %mean centering not needed here since you don't care about intercept
%    st_kernel = st_kernel(2:end); %remove intercept     
%     
%    %deconvolve
%    st_kernel = cat(1,zeros(win-1,1),-1*flipud(st_kernel));
%    st_kernel = st_kernel/max(st_kernel);
%    temp = deconvlucy(dff_probe(:,cur_probe),st_kernel);
% %    st_regressed = convn(padarray(st',[0,floor(length(st_kernel)/2)],'replicate','both')',flipud(st_kernel),'valid');           
%    
%     figure; plot(temp(1:1000)); yyaxis right;  hold on; plot(st(1:1000));
%     figure; plot(dff_probe(1:1000,cur_probe)); yyaxis right;  hold on; plot(st(1:1000));
% %%
%     figure; plot(xcorr(st-nanmean(st),temp-nanmean(temp),200,'normalized'))
%     hold on; plot(xcorr(st-nanmean(st),dff_probe(:,cur_probe)-nanmean(dff_probe(:,cur_probe)),200,'normalized'))
% %%
%    %predict st timing at timepoint t using following timepoints
%    win = round(opts.fps);
%    if mod(win,2)==0; win = win+1;  end 
%    predictors = createRollingWindow(dff_probe(2:end,cur_probe), win); %t-n:t-1
%    response = dff_probe(1:end-win,cur_probe);
%    auto_kernel = glmfit(predictors,response); %mean centering not needed here since you don't care about intercept
%    auto_kernel = auto_kernel(2:end); %remove intercept
%    %max normalize
% %    auto_kernel = auto_kernel/max(auto_kernel);
%    
%    %regress out autocorr
%    dff_regressed = convn(padarray(dff_probe(:,cur_probe)',[0,floor(length(auto_kernel)/2)],'replicate','both')',auto_kernel,'valid');           
%    figure; plot(dff_regressed);  hold on; plot(dff_probe(:,cur_probe));
%    %%
%    %predict imaging timing at timepoint t using preceding ephys
%    predictors = createRollingWindow(dff_probe(:,cur_probe), win);
%    response =  st(ceil(win/2):end-floor(win/2)); % get the middle timepoitn in window  
%    st_kernel = glmfit(predictors,response); %mean centering not needed here since you don't care about intercept
%    st_kernel = st_kernel(2:end); %remove intercept     
%     
%    %deconvolve   
%    st_kernel = st_kernel/max(st_kernel);
%    st_kernel = st_kernel-min(st_kernel);
%    temp = deconvlucy(dff_probe(:,cur_probe)-min(dff_probe(:,cur_probe)),st_kernel);
% %     temp = lucric(dff_probe(:,cur_probe)-min(dff_probe(:,cur_probe)),gp.d_gamma,gp.d_smooth,gp.d_kernel);
% %    temp = deconvreg(dff_probe(:,cur_probe),st_kernel);
% %    temp = convn(padarray(dff_probe(:,cur_probe)',[0,floor(length(st_kernel)/2)],'replicate','both')',st_kernel,'valid');           
%    
%     figure; plot(temp(1:1000)); yyaxis right;  hold on; plot(st(1:1000));
% %     figure; plot(dff_probe(1:1000,cur_probe)); yyaxis right;  hold on; plot(st(1:1000));
% figure; plot(xcorr(temp-nanmean(temp),st-nanmean(st),200,'normalized'));
% hold on; plot(xcorr(st-nanmean(st),st-nanmean(st),200,'normalized'));
% %%   %fit your kernel | predict imaging timing at timepoint t using preceding ephys  
%    win = 20;%round(opts.fps);
%    if mod(win,2)==0; win = win+1;  end 
%    predictors = createRollingWindow(dff_probe(:,cur_probe), win);
%    response =  st(ceil(win/2):end-floor(win/2)); % get the middle timepoitn in window  
%    kernel = glmfit(predictors,response); %mean centering not needed here since you don't care about intercept
%    kernel = kernel(2:end); %remove intercept     
%    kernel = kernel/max(kernel); %normalize
%      
%    %estimate spiking using foopsi
%    [c, s, options] = deconvolveCa(dff_probe(:,cur_probe), 'kernel',-1*flipud(kernel),'foopsi');
%    figure; hold on; 
%    plot(dff_probe(1:1000,cur_probe),'color','k');
%    plot(c(1:1000),'color','r');
%    plot(s(1:1000),'color','b');
%    plot(st(1:1000),'color','c');
%    figure; hold on; plot(xcorr(st-nanmean(st),s-nanmean(s),200))
% 
% 
%     %best using the LR
%    st_kernel = kernel;   
%    st_kernel = st_kernel/max(st_kernel);
%    st_kernel = st_kernel-min(st_kernel);
%    temp = deconvlucy(dff_probe(:,cur_probe)-min(dff_probe(:,cur_probe)),st_kernel);
%    
%    figure; hold on; 
%    plot(dff_probe(1:300,cur_probe),'color','k');
%    plot(st(1:300),'color','c');
%   yyaxis right
%    plot(temp(1:300),'color','r');   
%    legend('dff','predicted','spiking')
%    title('LR-deconv example trace')
%    
%    figure; hold on; plot(xcorr(dff_probe(:,cur_probe),dff_probe(:,cur_probe),200,'normalized'),'color','k')
%    plot(xcorr(nn_output-nanmean(nn_output)',200,'normalized'),'color','r')
%    plot(xcorr(st-nanmean(st),200,'normalized'),'c')
%    plot(xcorr(st-nanmean(st),temp-nanmean(temp)',200,'normalized'),'color','g')
%    plot(xcorr(st-nanmean(st),dff_probe(:,cur_probe)-nanmean(dff_probe(:,cur_probe)),200,'normalized'),'color','m')
%    legend('dff-dff','pred-pred','spike-spike','pred-spike','dff-spike')
%    title('LR-deconv auto/cross corr')   
%    figure; hold on; 
%    subplot(3,1,1); histogram(st); title('true fr'); 
%    subplot(3,1,2); histogram(temp); title('predicted fr')
%    subplot(3,1,3); histogram(dff_probe(:,cur_probe)); title('dff')
%    title('FR distribution')   
%    
%    %get the kernel
%    st_regressed = convn(padarray(st',[0,floor(length(st_kernel)/2)],'replicate','both')',st_kernel,'valid');           
%   
%    %img --> ap kernel | predict fr at each timepoint using 1s imaging data    
%    %individual points
%    win = opts.fps*2;
%    if mod(win,2)==0; win = win+1;  end %odd so can appropriately center response within window
%    predictors = createRollingWindow(dff_probe(:,cur_probe), win);
%    response =  st(ceil(win/2):end-floor(win/2)); % get the middle timepoitn in window  
% %    predictors = createRollingWindow(st, win);
% %    response =  dff_probe(ceil(win/2):end-floor(win/2),cur_probe); % get the middle timepoitn in window     
%    kernel = NaN(numel(win),1);
%    for i = 1:size(predictors,2)
%       temp = glmfit(predictors(:,i),response); %mean centering not needed here since you don't care about intercept
%       kernel(i,1) = temp(2);
%    end
% %    kernel = kernel(2:end); %remove intercept. 
%       
%    figure; hold on;
%    t=[-1*floor(win/2):floor(win/2)]/30;   
%    plot(t,kernel,'color','k','linewidth',2);   
%    yval = get(gca,'ylim');   
%    line([0 0],[yval],'color','k','linestyle',':','linewidth',2)  
%    xlim([-0.75 0.75]); ylim([yval]);
%    title(sprintf('Probe %d Kernel Single Points',cur_probe));
%    xlabel('flourescence Lag (s)');    
%    
%       %kernel full window
%    if mod(win,2)==0; win = win+1;  end %odd so can appropriately center response within window
%    predictors = createRollingWindow(dff_probe(:,cur_probe), win);
%    response =  st(ceil(win/2):end-floor(win/2)); % get the middle timepoitn in window  
%    kernel= glmfit(predictors,response); %mean centering not needed here since you don't care about intercept
%    kernel = kernel(2:end); %remove intercept. 
%       
%    figure; hold on;
%    t=[-1*floor(win/2):floor(win/2)]/30;   
%    plot(t,kernel,'color','k','linewidth',2);   
%    yval = get(gca,'ylim');   
%    line([0 0],[yval],'color','k','linestyle',':','linewidth',2)  
%    xlim([-0.75 0.75]); ylim([yval]);
%    title(sprintf('Probe %d Kernel Full Window',cur_probe));
%    xlabel('flourescence Lag (s)'); 
%     
%    %test how well it works by looks at the autocorrelation before and after
%    dff_conv = convn(padarray(dff_probe(:,cur_probe)',[0,floor(length(kernel)/2)],'replicate','both')',kernel,'valid');
%    figure; hold on; 
%    subplot(2,3,1)
%    histogram(dff_probe(:,cur_probe));
%    subplot(2,3,2)
%    histogram(dff_conv(:));
%    subplot(2,3,3)
%    histogram(st(:));
%    subplot(2,3,[4,6]); hold on;
%    plot(dff_probe(1:opts.fps*30,cur_probe),'color',[0.5 0.5 0.5 0.5],'linewidth',2)
%    plot(dff_conv(1:opts.fps*30),'color','k','linewidth',2)
%    yyaxis right
%    plot(st(1:opts.fps*30),'color','r','linewidth',2,'linestyle',':')
%    
%    %deconvolved img --> ap xcorr and deconv img --> deconv --img
%    [xcorrkernel_deconv,t] = xcorr(st-nanmean(st),dff_conv'-nanmean(dff_conv),maxlag,'normalized');   
%    xcorrkernel = xcorr(st-nanmean(st),dff_probe(:,cur_probe)'-nanmean(dff_probe(:,cur_probe)),maxlag,'normalized');
%    deconv_auto = xcorr(dff_conv'-nanmean(dff_conv),maxlag,'normalized');
%    img_auto = xcorr(dff_probe(:,cur_probe)'-nanmean(dff_probe(:,cur_probe)),maxlag,'normalized');
%    figure; hold on; 
%    p0=plot(t/30,xcorrkernel_deconv,'k','linewidth',2);    
%    p1=plot(t/30,xcorrkernel,'k','linewidth',2,'linestyle',':');    
%    xval = get(gca,'xlim');
%    yval = get(gca,'ylim');
%    line([0 0],[-.1 .1],'color','k','linestyle',':')   
%    xlabel(' Lag (s)')
%    ylabel('Correlation');
%    xlim([xval]); ylim([yval]);   
%    yyaxis right
%    p2=plot(t/30,deconv_auto,'g','linewidth',2,'linestyle','--');
%    p3=plot(t/30,st_auto,'r','linewidth',2,'linestyle',':');
%    p4=plot(t/30,img_auto,'b','linewidth',2,'linestyle','--');
%    legend([p0,p1,p2,p3,p4],'deconv-->st','img-->st','deconv-->deconv','st-->st','img-->img');
%    set(gca,'ycolor','k'); ylabel('Autocorrelation');
%    
%    %include in title the quality of deconv fit to fr
%    title(sprintf('Probe %d Xcorr Deconv RMSE=%0.2g Img RMSE=%0.2g',cur_probe,sqrt(mse(st-dff_conv)),sqrt(mse(st-dff_probe(:,cur_probe)))));
%       
%    %ap --> ap autokernel img --> img autokernel 
%    predictors = createRollingWindow(st, win);
% %    predictors(:,ceil(win/2)-5:ceil(win/2)+5)=0;
%    response =  st(ceil(win/2):end-floor(win/2)); % get the middle timepoitn in window   
%    kernel = glmfit(predictors,response);
%    kernel = kernel(2:end); %remove intercept. 
%       
%    figure; hold on;
%    t=[-1*floor(win/2):floor(win/2)]/30;   
%    plot(t,kernel,'color','k','linewidth',2);   
%    yval = get(gca,'ylim');   
%    line([0 0],[yval],'color','k','linestyle',':','linewidth',2)  
%    xlim([-0.75 0.75]); ylim([yval]);
%    title(sprintf('ST GLM d%d r%d MUA%d',params.depth,params.radius,params.mua));
%    xlabel('flourescence Lag (s)');      
%    
% end %probe loop
% 
% %save off figures
% handles = get(groot, 'Children'); 
% file_list={};
% for i = 1:numel(handles)
%     file_list{i} = [pwd filesep ,sprintf('%d_kernel.pdf',i)];
%     print(handles(i),file_list{i},'-dpdf','-bestfit');
% end
% append_pdfs(sprintf('pca_kernels_muas%d_d%d_r%d.pdf',params.mua,params.depth,params.radius),file_list{:});
% 
% % %delete original figures
% cellfun(@(x) delete(x),file_list)
% 
% 
% end %function
% 
% % add in self prediction for both imaging and neural
% % remake so using the entire data
% % add probe adjustments
% 
% 
% %%old
% % num_iter =1000;  
% % binsize = opts.fps*2.5;
% % xcorrkernel = NaN(num_iter,binsize*2+1);   
% % dur=(opts.fps*10*60);       
% % for i = 1:num_iter   
% %   idx = randi([dur+1,numel(st)-dur]);
% %   [xcorrkernel(i,:),t] = xcorr(st(idx:idx+dur),dff_probe(idx:idx+dur,cur_probe)',binsize,'normalized');
% % end
% 
% 
% %    binsize = opts.fps*2;
% %    kernel=[];
% %    for j = 1:1000
% %        str = randi([2*binsize,numel(st)-(2*opts.fps*60)]);
% %        count=1;
% %        predictors = [];
% %        response = [];   
% %        for i = str:str+opts.fps*60
% %            predictors(count,:)= dff_probe(i-binsize:i+binsize,cur_probe)';
% %            response(count) = st(i);    
% %            count= count+1;
% %        end
% %        kernel(j,:) = glmfit(predictors,response);
% %    end
% %    %remove intercept
% %    kernel = kernel(:,2:end);
% %    %maximum normalized
% %    kernel_norm = kernel./max(abs(kernel),[],2);
% %    kernel_mean = nanmean(kernel_norm,1);       
% %    kernel_mean = kernel_mean./sum(kernel_mean.^2);