function CompareSubspaceUniformity(randseed,verbose)
%Camden MacDowell - timeless
%get the unifromity (avg entropy) of contributions to a brain area from
%other areas. 
%compares to a shuffled distribution

if nargin <2; verbose = 0; end
%Add paths
if ispc
    addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\fpCNMF'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Widefield_Imaging_Analysis'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\GenUtils'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Ephys'));
    % Save off values
    savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\SubspaceUniformity';
else
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/fpCNMF'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/GenUtils'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Ephys'));
    % Save off values
    savedir = '/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/Analysis/SubspaceUniformity';
end

data = LoadSubspaceData('in_grouped');
rng(randseed);
if randseed <1
    fprintf('\nShuffle time!');
end
[~, area_all] = LoadVariable(data,'rrr_beta','VIS',1); %load areas
get_ent = @(p) -sum(p.*log2(p));
sweepidx = 0:0.25:1;
ent = NaN(10,8,84,numel(sweepidx));
strongest = NaN(10,8,84,numel(sweepidx),8);
tic
for cur_d = 1:10
    fprintf('\nworking on dim %d of %d',cur_d,10);
    for cur_area = 1:numel(area_all)
        targ_area= area_all{cur_area};       
        [ball, ~] = LoadVariable(data,'rrr_beta',targ_area,cur_d); %load betas
        area=LoadVariable(data,'beta_region',targ_area);                
        ball = reshape(ball,size(ball,1)*size(ball,2),size(ball,3));
        area = reshape(area,size(area,1)*size(area,2),size(area,3));

        %loop through all models
        for cur_model = 1:size(ball,1)
            b = abs(ball(cur_model,:));
            a = area(cur_model,:);
            [~,idx] = sort(b,'descend');
            b = b(idx);
            a = a(idx);
            %remove nan
            badidx = isnan(b);
            a(badidx)=[];
            b(badidx)=[];
            if randseed >1
                a = a(randperm(numel(a),numel(a))); %randomly permute the area labels
            end            

            %sweep fractions            
            for i = 1:numel(sweepidx)
                areatemp = a(1:ceil(sweepidx(i)*numel(b)));
                y = arrayfun(@(n) nansum(areatemp==n),unique(areatemp),'UniformOutput',1); %get the numbers for each area
                %convert to fraction of that areas total number of neurons
                yall = arrayfun(@(n) nansum(a==n),unique(areatemp),'UniformOutput',1);
                y = y./yall;     
                %get the entropy of that distribution
                ent(cur_d,cur_area,cur_model,i) = get_ent(y);
                %get the strongest contributor
%                 temp = unique(areatemp);
                strongest(cur_d,cur_area,cur_model,i,unique(areatemp)) = y;
            end
        end
    end
end

save([savedir, filesep, sprintf('AbsBetas_run%d.mat',randseed)],'strongest','ent');

% Plot the strongest connections
if verbose == 1
    fp = fig_params_cortdynamics;
    area_label = area_all;
    x = permute(strongest,[2,3,1,5,4]);
%     x = x(:,:,1,:,:)
    x = reshape(x,size(x,1),size(x,2)*size(x,3),size(x,4),size(x,5));
    y = NaN(size(x,1),size(x,2),size(x,3));
    for i = 1:size(x,1)
        for j = 1:size(x,2)
            for z = 1:size(x,3)
                temp = [0,squeeze(x(i,j,z,2:end-1))',1];
    %             y(i,j,z) = nanmean(squeeze(x(i,j,z,2:end-1)));
                y(i,j,z) = trapz(temp/numel(temp));
            end
        end
    end

    w = [1,0.6 0.25];
    yy = squeeze(nanmean(y,2));
    %get the top contributions to each
    for i = 1:size(yy,1)
       temp = yy(i,:);
       [~,idx] = maxk(temp,numel(w));
       for j = 1:numel(w)
           temp(idx(j))=w(j);
       end
       temp(~ismember(1:numel(temp),idx))=0;
       yy(i,:) = temp;
    end
    figure; hold on;
    circularGraph(yy,'colormap',fp.c_area,'label',area_label,'normVal',0.4); %plot
    set(findall(gca, 'type', 'text'), 'visible', 'on') 
    fp.FigureSizing(gcf,[3 3 6 6],[10 10 14 10])

    for i = 1:size(yy,1)
        a = yy;   
        a(~ismember(1:size(a,1),i),:)=0;
        figure; hold on;
        circularGraph(a,'colormap',fp.c_area,'label',area_label,'normVal',0.4); %plot
        set(findall(gca, 'type', 'text'), 'visible', 'on') 
        fp.FigureSizing(gcf,[3 3 6 6],[10 10 14 10])
    end
    
    %% Plot the top relationship
    w = [1];
    yy = squeeze(nanmean(y,2));
    %get the top contributions to each
    for i = 1:size(yy,1)
       temp = yy(i,:);
       [~,idx] = maxk(temp,numel(w));
       for j = 1:numel(w)
           temp(idx(j))=w(j);
       end
       temp(~ismember(1:numel(temp),idx))=0;
       yy(i,:) = temp;
    end
    for i = 1:size(yy,1)
        a=yy;
        a(~ismember(1:size(a,1),i),:)=0;
        figure; hold on;    
        circularGraph(a,'colormap',fp.c_area,'label',area_label,'normVal',0.4); %plot
        set(findall(gca, 'type', 'text'), 'visible', 'on') 
        fp.FigureSizing(gcf,[3 3 6 6],[10 10 14 10])
        title('strongest');
    end
    
    %plot the second relationship
    w = [1,0.7];
    yy = squeeze(nanmean(y,2));
    %get the top contributions to each
    for i = 1:size(yy,1)
       temp = yy(i,:);
       [~,idx] = maxk(temp,numel(w));
       idx = idx(2);       
       temp(idx)=w(2);
       temp(~ismember(1:numel(temp),idx))=0;
       yy(i,:) = temp;
    end
    for i = 1:size(yy,1)
        a=yy;
        a(~ismember(1:size(a,1),i),:)=0;
        figure; hold on;    
        circularGraph(a,'colormap',fp.c_area,'label',area_label,'normVal',0.4); %plot
        set(findall(gca, 'type', 'text'), 'visible', 'on') 
        fp.FigureSizing(gcf,[3 3 6 6],[10 10 14 10])
        title('second strongest');
    end
    
    %plot the thrid relationship
    w = [1,0.7,0.5];
    yy = squeeze(nanmean(y,2));
    %get the top contributions to each
    for i = 1:size(yy,1)
       temp = yy(i,:);
       [~,idx] = maxk(temp,numel(w));
       idx = idx(3);       
       temp(idx)=w(3);
       temp(~ismember(1:numel(temp),idx))=0;
       yy(i,:) = temp;
    end
    for i = 1:size(yy,1)
        a=yy;
        a(~ismember(1:size(a,1),i),:)=0;
        figure; hold on;    
        circularGraph(a,'colormap',fp.c_area,'label',area_label,'normVal',0.4); %plot
        set(findall(gca, 'type', 'text'), 'visible', 'on') 
        fp.FigureSizing(gcf,[3 3 6 6],[10 10 14 10])
        title('third strongest');
    end    
    
    %% Add plots
    fp = fig_params_cortdynamics;
    %plot entopy across dimensions
    % x = nanmean(ent,[2:4]);
    %replace zeros with NaN (since missing areas)
    ent(ent==0)=NaN;
    x = permute(ent,[2,3,1,4]);
    % x = squeeze(x(:,:,:,:));
    x = x(:,:,:,2:end-1);
    %average across dimensions
    % x = squeeze(nanmean(x,[3]))
    x = reshape(x,size(x,1),size(x,2)*size(x,3)*size(x,4));
    x = pairedBootstrap(x',@nanmean); % bootstrap the distributions
    % x = x';
    %reorder 
    [~,idx] = sort(nanmean(x),'descend');
    x = x(:,idx);   

    % plot as a violin showing each region
    col = fp.c_area;
    col = arrayfun(@(n) col(n,:),1:size(col,1),'UniformOutput',0);
    col = col(idx);
    area_name = area_all(idx);

    %flatten per area
    figure; hold on;
    vp = CompareViolins(x',fp,'col',col,'connectline',[],'plotspread',0,'divfactor',0.5);
    fp.FigureSizing(gcf,[3 2 6 4],[10 10 20 10]); 
    set(gca,'XTickLabel',area_name,'XTickLabelRotation',45)
    fp.FormatAxes(gca); box on; grid on
    ylabel('Sub D / Local D');
    
    savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\SubspaceUniformity';
    saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'SubspaceUniformityPlots',savedir,0); close all
end


%% If you want to evaluate the final data
%load the entropy and average across all dimensions and sampled points









end %function end









