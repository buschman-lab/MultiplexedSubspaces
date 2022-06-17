function [rho, area_label, pval] = CompileInteractions(data,useProj,ndim,type)
%Camden MacDowell - timeless
%NOTE; beta direction is random, so this returns absolute value of
%correlations

if nargin <2; useProj = 1; end %def = 1; use the timepoint/trial projection along each beta (slower). 0 = just correlate betas
if nargin <3; ndim = 4; end %number of dimensions to use
if nargin <4; type = 'uniform'; end % do an exhaustive search or just get the same dimension in all areas. 

nrec = numel(data); %num  recordings
nmotif = size(data{1},2); %num motifs;

%full list of areas (recs 1,2,5,6 have all areas; 3 is missing RSP and BFD. 4 if missing RSP)
area_label = data{1}(1).area_label; 

rho = cell(nrec,nmotif); %using cell arrays to help keep things straight since lots of same sizes (8 brain regions)
pval = cell(nrec,nmotif); %using cell arrays to help keep things straight since lots of same sizes (8 brain regions)
for cur_rec = 1:nrec %loop through recordings
    fprintf('\nworking on rec %d',cur_rec);
    for cur_motif = 1:nmotif %loop through motifs    
        %load betas
         betas = loadBetas(data,cur_rec,cur_motif,ndim,area_label,useProj);         
         %get correlations
         [rho{cur_rec,cur_motif}, pval{cur_rec, cur_motif}] = getRho(betas,type);         
    end %motif loop    
end %rec loop 

end %function 


function [rho, pval] = getRho(betas,type)
    %returns rho, which is a 8x1 cell array for all brain regions. 
    %within each cell is a 8x8xndim similarity matrix of how similar
    %regions are represented in a brain region
    rho = cell(size(betas,1),1); 
    pval = cell(size(betas,1),1); 
    for cur_a = 1:size(betas,1) %loop through each row of betas (predictor area - see loadBetas)
       b = betas(cur_a,:);
       %populate with NaN of appropriate size
       idx = find(cellfun(@(x) ~isempty(x),b),1,'first');
       if ~isempty(idx) %missing brain region in this rec
           ndim = size(b{idx},2);
           nneu = size(b{idx},1);
           b(cellfun(@isempty,b))={NaN(nneu,ndim)};

           switch type
               case 'uniform' %same dimensions in all areas
                   %loop through each target dimension 
                   r = NaN(8,8,ndim);
                   for cur_d = 1:ndim                    
                       bb = cellfun(@(x) x(:,cur_d),b,'UniformOutput',0);
                       bb = cat(2,bb{:});
                       a = abs(corr(bb,'type','pearson')); %this is the correlation between projections/betas between each target region in area cur_a for dimension cur_d
                       a(1:1+size(a,1):end) = NaN;%set diag to NaN      
                       r(:,:,cur_d) = a; 
                   end
               case {'exhaustive', 'optimal'} %correlation between all dimensions in all regions
                   %the 4 dimensions of the target will be all NaN (since it can't represent itself)                 
                   bb = cat(2,b{:}); 
                   [cur_rho, cur_p] = corr(bb, 'type', 'pearson');
                   a = abs(cur_rho); 
                   a(1:1+size(a,1):end) = NaN; %the 4 dimensions of the target will be all NaN   
                   r = a; 
                   p = cur_p;
                   %remove correlations between dimensions of the same target regions
%                    x = 1:ndim:size(r,1);
%                    for i = 1:numel(x)
%                        r(x(i):(x(i)+ndim-1),x(i):(x(i)+ndim-1))=NaN;
%                    end
               otherwise 
                   error('unknown correlation approach');
           end %switch 
           rho{cur_a} = r; 
           pval{cur_a} = p;
       end %empty if
    end %cur_a loop
    
    %fill any missing brain areas with NaN to allow easy averaging
    padsize = size(rho{find(cellfun(@(x) ~isempty(x),rho),1,'first')});
    rho(cellfun(@isempty,rho))={NaN(padsize)};
    pval(cellfun(@isempty,pval))={NaN(padsize)};
    
end


function [betas,targidx] = loadBetas(data,cur_rec,cur_motif,ndim,area_label,useProj)
    %retruns betas which is an 8x8 cell array. Each cell cointains matrix
    %of betas weights of a region (nuerons x dimensions)
    % useProj, then returns (timepoints x trials) x dimensions matrix of
    % activtiy in a region projected along those beta weights
    %rows of beta are predictor regions. columns are targets
    %i.e., row 1 of beta column 2 is (neurons x dimensions) betas for Hipp predicting area MOs
    %Gutchecks: 
    %betas should have nothing on the diagnol
    %betas for rec three rows and columns 4/6 should be empty and rec 4 row/column 3

    b = data{cur_rec}(cur_motif).rrr_B; %each cell is beta weights of n-1 regions predicting each target region
    cur_areas = data{cur_rec}(cur_motif).area_label; %the list of target regions
    
    %just keep desired dimensions
    b = cellfun(@(x) x(:,1:ndim),b,'UniformOutput',0);     
    
    %reorganize betas so each cell is the betas for a given area when predicting all other areas
    g = data{cur_rec}(cur_motif).grouping; %each cell is region of beta weights of n-1 regions predicting each target region             
    betas = cell(numel(area_label),numel(area_label));
    targidx = ismember(area_label,cur_areas); %get true index for targets (out of 8 areas, to match across recs)
    for cur_a = 1:numel(cur_areas) 
        idx = strcmp(area_label,cur_areas(cur_a)); %get true index (out of 8 areas, to match across recs)        
        if ~isempty(idx)
            betas(idx,targidx) = arrayfun(@(n) b{n}(g{n}==cur_a,:),1:numel(g),'UniformOutput',0);  
        end
    end  
    
    %if projecting 
    if useProj == 1
       y = data{cur_rec}(cur_motif).area_val;
       y = cellfun(@(x) normalizeToBaseline(x,[1:2],'mean'), y,'UniformOutput',0);%normalize to baseline
       y = cellfun(@(x) x(:,3:end,:), y,'UniformOutput',0);%use non-baseline
       y = cellfun(@(x) x-nanmean(x,3), y,'UniformOutput',0);%subtract the psth
       y = cellfun(@(x) reshape(x,[size(x,1),size(x,2)*size(x,3)])', y,'UniformOutput',0);%reshape into trial/timepoint x neuron   
       
       %reshape to match predictors in beta
       yy = cell(8,1);
       yy(targidx)=y;
       yy = repmat(yy,1,8);
       %temporily replace emptys with nan for multiplication
       badidx = cellfun(@isempty,betas);
       yy(badidx) = {NaN};
       betas(badidx) = {NaN};
       betas = cellfun(@(x,y) x*y, yy,betas,'UniformOutput',0);
       betas(badidx) = {[]};
       
    end %proj if

end




