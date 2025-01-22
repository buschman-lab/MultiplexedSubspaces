%Figure generation code for 
%Multiplexed Subspaces Route Neural Activity Across Brain-wide NetworksMultiplexed Subspaces Route Neural Activity Across Brain-wide Networks
%MacDowell et al., 
%Corresponding author: tbuschma@princeton.edu

load('area_label.mat');

%% Figure 1;

%1F
load('Figure1f.mat');
Figure1f(mdl_perf)

%1H
load('Figure1h.mat');
Figure1h(rrr_mdl,area_label)

%1I
load('Figure1i.mat');
Figure1i(Dshared,Dlocal)

%1J
load('Figure1J.mat')
Figure1j(x,col,area_label) 


%% Figure 2;
 %2C
load('Figure2C.mat')
Figure2c(b,area,area_label); 

%2D
load('Figure2D.mat')
Figure2d(x,ci,area_all)

%2F
load('Figure2F.mat')
Figure2f(behav_rho)

%% Figure 4

%4C
load('Figure4C.mat')
Figure4c(ovrlap)

%% Figure 5

%5a
load('Figure5A.mat');
Figure5(clust_colmap,clust_size,k,sorted_sim_mat)

%% Figure 6

%6A 
load('Figure6a.mat')
Figure6a(y) 

%6B (example hexagon plot creation)
load('Figure6B.mat');
Figure6b(rho,sig_thresh)

%% Figure 7

%7A
load('Figure7a.mat')
Figure7a(a,b) 

%7B
load('Figure7b.mat')
Figure7b(a,b,area)  

%% Figure 8

%8A
load('Figure8a.mat')
Figure8(theta,thetacon)

%8B
load('Figure8b.mat')
Figure8b(theta,ci)

%8C
load('Figure8c.mat')
Figure8c(theta,xpos)

%% Figure 9

%9A
load('Figure9a.mat')
Figure9a(proj1,proj2)
















