function [dff_z, cursor, target_hit, c1_bool, c2_val, c2_bool, c3_val, c3_bool] = ...
    dff2cursor_target(dff, bData, cursor_zscore_bool) 
%4.19.19
%INPUT: smoothed dff
%1) z-score dff
%2) calculate cursor value
%3) determine if target is hit
% OUTPUT: 
% c1: cursor
% c2: E1_mean > E1_thresh
% c3: E2_subord_mean > E2_subord_thresh
%
% relevant fields from bData: 
%n_mean, n_std, decoder, E1_sel_idxs, E2_sel_idxs, E1_thresh, E2_subord_thresh, T)
%decoder: num_neurons x 1

%E2_subord_thresh: size (num_E2-1 x 1)
%cond1: neural*decoder >= T
%cond2: E1<= E1_thresh
%cond3: E2_subord_mean >= E2_subord_thresh(E2_dom)

%neural: 1 x num_neurons
%decoder: num_neurons x 1
%E_id: num_neurons x 1
%E1_thresh: 1x1
%E2_subord_thresh: 1xnum_E2
%T: 1x1

%Ensemble Identity Info:
% E1_sel = E_id==1; 
% E1_sel_idxs = find(E1_sel); 
% num_E1 = length(E1_sel_idxs); 
% 
% E2_sel = E_id==2; 
% E2_sel_idxs = find(E2_sel); 
% num_E2 = length(E2_sel_idxs); 

num_E2 = length(bData.E2_sel_idxs); 

%z-score:
dff = dff(:).';
dff_z = dff; 
if(cursor_zscore_bool)
    dff_z = (dff-bData.n_mean)./bData.n_std;
    dff_z = dff_z(:).'; %set dff_z to be a row    
    n_analyze = dff_z;
else
    n_analyze = dff;
end

E1 = n_analyze(bData.E1_sel_idxs); 
E2 = n_analyze(bData.E2_sel_idxs); 

%c1: cursor
cursor = n_analyze*bData.decoder;
c1_val = cursor; 

%c2: E1_mean
E1_mean = mean(E1); 
c2_val = E1_mean;

%c3: E2_subord
E2_sum                          = sum(E2); 
[E2_dom_samples, E2_dom_sel]    = max(E2, [], 2); %E2_dom  CAREFUL THIS MAY BRING TWO!!
% E2_dom_samples = n_analyze(E2_dom_sel(1));
E2_subord_mean                  = (E2_sum - E2_dom_samples)/(num_E2-1);

c3_val  = E2_subord_mean;

c1_bool = cursor >= bData.T1; 
c2_bool = E1_mean <= bData.E1_thresh;
c3_bool = E2_subord_mean >= bData.E2_subord_thresh(E2_dom_sel);

%target hit
% target_hit = c1_bool & c2_bool & c3_bool; 
target_hit = c1_bool;
