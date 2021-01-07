% std_envtopo() - Creates an envtopo() plot for STUDY component clusters. Each cluster ERP
%                 is normalized by the number of ICs per group. When generating multiple 
%                 plots for each combination of conditions, cluster line colors
%                 are shared across plots.
% 
% statPvaf()    - subfunction for bootstrap statistics for p-values and confidence interval.
%                 stat_surrogate_pvals and stat_surrogate_ci are used.
%
% Usage:
%           >> std_envtopo(STUDY, ALLEEG, 'key1', 'val1', ...);
%
%           When using with GUI, choose time range (in ms) to plot, time range
%           (in ms) to rank cluster contributions, and specify cluster(s) of interest
%           or number of largest contributing cluster(s) to plot. envtopo() plots can 
%           be generated for a single variable combination (using "STUDY variable
%           combination"), a subtraction between two variable combinations (using
%           "Subtraciton"), or an interaction between four variable combinations 
%           (using "Interaction"). statPvaf is accessed by entering
%           [alpha numIterations] into "pvaf statistics" field in GUI. 
%
% Inputs:
%
%   STUDY   an EEGLAB STUDY structure containing EEG structures
%
%   ALLEEG  the ALLEEG data structure; can also be an EEG dataset structure.
%
% Optional inputs:
%
%  'baseline' [integer ms integer ms] - a new baseline (beginning and ending, respectively)
%             to remove from the grand and cluster ERPs.
%
%  'clustnums' [integer array] vector of cluster numbers to plot.  Else if
%              int < 0, the number of largest contributing clusters to plot
%              {default|[] -> 7}
%          
%  'downsampled' ['on'|'off'] Topography downsampled from 67x67 to 23x23
%                Always 'on' for statPvaf calculations. {Default 'off'} 
%
%  'drawvertline' [integer ms] vector of times at which to plot drawvertlineical dashed lines
%                   {default|[] -> none}
%
%  'fillclust' [integer] fill the numbered cluster envelope with red. {default|[]|0 -> no fill}
%
%  'fillclustcolor' [float float float] where inputs are 0<=value<=1. to create a color
%                    to fill the summed selected clusters. {dafault [1 0 0]}
%
%  'fillenvsum' ['selective'|'all'|'off'] fill or not the envelope of the
%               summed selected cluster. 'select' fills only limitcontribtime
%               time range. 'all' fills all timerange. 'off' does not fill.
%               See also 'fillenvsumcolor' so choose a filling color. {default: 'selective'}
%
%  'fillenvsumcolor' [float float float] where inputs are 0<=value<=1. To create a color
%                      to fill the summed selected clusters. {dafault [0.875 0.875 0.875]}
%
%  'forcePositiveTopo' ['on'|'off'] flip scalp topo polarities to negative if
%                      abs(min(envelope)) > abs(max(envelope)). {default:'off'}
%
%  'limitamp' [float uV float uV]. Specify minimum and maximum amplitude limit respectively
%              {default: use data uV limits}
%
%  'limitcontribtime' [integer ms integer ms]  time range in which to rank cluster contributions
%                      (boundaries = thin dotted lines) {default|[]|[0 0] -> plotting limits}
%  
%  'lpf' [float Hz] low-pass filter the ERP data using the given upper pass-band edge
%                   frequency.  Filtering by pop_eegfiltnew(). {default: no lfp}
%  
%  'normalization'   ['normIC'|'normSubj'|'off'] Normalize projection by number of ICs, number of
%                    subjects, or no normalization, respectively. {Default 'normIC'} 
%
%  'onlyclus' ['on'|'off'] dataset components to include in the grand ERP.
%             'on' will include only the components that were part of the
%              clustering. For example, if components were rejected from
%              clustering because of high dipole model residual variance,
%              don't include their data in the grand ERP.
%             'off' will include all components in the datasets except
%              those in the subtructed ('subclus') clusters {default 'on'}.
%
%  'measure' ['pvaf'(default)|'ppaf'|'maxp'|'rltp'|'auc']
%            When called from GUI, only PVAF and PPAF are supported.
%            'pvaf', sort components by percent variance accounted for (eeg_pvaf())
%            pvaf(comp) = 100-100*mean(var(data - back_proj))/mean(var(data));
%            'ppaf', sort components by percent power accounted for (ppaf) 
%            ppaf(comp) = 100-100*Mean((data - back_proj).^2)/Mean(data.^2);
%            %%%%%%%%%%%%%%%% Not suppported by GUI %%%%%%%%%%%%%%%%%%%%
%            'maxp', maximum power
%            maxp(comp) = max(sum(cluster(:,selectedtimewindow)^2));
%            'rltp', sort components by relative power 
%            rltp(comp) = 100*Mean(back_proj.^2)/Mean(data.^2);
%            'auc', sort components by area under a curve (AUC) ratio.
%            auc(comp) = 100*(AUC of envelope by a selected cluster)/(auc of outermost envelope)
%
%  'normalizeClusterPvaf' ['on'|'off'] normalizes clusterPvaf so that the sum of clusterPvaf is 100%. {Default 'off'} 
%
%  'subclus' [integer array] vector of cluster numbers to omit when computing
%            the ERP envelope of the data (e.g., artifact clusters).
%            By default, these clusters are excluded from the grandERP envelope.
%            See the option 'onlyclus'.
%
%  'timerange' [integer ms integer ms] data epoch start and end input
%              latencies respectively. {default: epoch length}
%
%  'topotype' ['inst'|'ave'] If 'inst', show the instantaneous map at the specific timepoint
%             specified by the line. If 'ave', show that of the mean which is the same as the
%             stored topomap. {default:'inst'}
%
% See also: eegplugin_std_envtopo, pop_std_envtopo, std_envtopo, envtopo, stat_surrogate_pvals, stat_surrogate_ci
% 			Lee et al. "Non-parametric group-level statistics for source-resolved ERP analysis." Engineering 
%				in Medicine and Biology Society (EMBC), 2015 37th Annual International Conference of the IEEE. IEEE, 2015.

% Author: Makoto Miyakoshi, 2011- JSPS;INC,SCCN,UCSD.
%         Former version of this function was developed by Hilit Serby, Arnold Delorme, and Scott Makeig.
% History:
% 01/07/2020 Makoto. Stats revisited. Now only one cluster selection is allowed to simplify the stats.
% 06/27/2019 Makoto. Improved efficiency.
% 09/19/2018 Makoto. Conversion from numeric condition labels to strings supported. Converted from double to single. 'Else plot these cluster numbers only' debugged. fillclust option fixed.
% 08/14/2018 Makoto. Reorganized envtopo_plot(). PVAF/PPAF switch implemented.
% 08/10/2018 Makoto. Updated. Compatible with EEGLAB15.
% 04/06/2016 ver 3.00 by Clement.statPvaf released. std_envtopo updated with GUIDE based GUI. Option 'normalization' and 'downsampled' added; 'diff' removed (obsolete because of GUI). Help updated.
% 01/19/2016 ver 2.42 by Makoto. Changed 1:size(setinds{1,1},2) to 1:size(STUDY.group,2) Thanks Alvin Li!
% 01/08/2016 ver 2.41 by Makoto. Removed 'Store the precomputed results' removed.
% 11/25/2015 ver 2.40 by Makoto. Session bug fixed.
% 02/09/2015 ver 2.39 by Makoto. Subtraction figure bug fixed.
% 01/30/2015 ver 2.38 by Makoto. Figure overlap problem fixed.
% 12/26/2014 ver 2.37 by Makoto. Bug fixed (when no group, plotting failed). Axis labels renewed. Low-pass filter order shortened to be a half.
% 08/04/2014 ver 2.36 by Makoto. Reports the number of subjects and ICs for each group.
% 07/08/2014 ver 2.35 by Makoto. Outlier cluster detection updated (now with try...catch...)
% 03/21/2014 ver 2.34 by Makoto. Help updated.
% 02/13/2014 ver 2.33 by Makoto. Debug on fillclust. Added fillclustcolor. Fine-tuned plot layout and design.
% 02/06/2014 ver 2.32 by Makoto. Commneted out the STUDY update at the end of the main function.
% 02/06/2014 ver 2.31 by Makoto. Bug fix in the eegplugin_std_evntopo.
% 02/04/2014 ver 2.30 by Makoto. Pvaf calculation now uses mean(var(ERP... for compatibility with envtopo().
% 02/03/2014 ver 2.20 by Makoto. Session supported (can't be used with group)
% 11/14/2013 ver 2.14 by Makoto. Bug if no group fixed (line 340). Use of lpf option forces to recompute data.
% 10/11/2013 ver 2.13 by Makoto. pvaf: xx colon removed.
% 10/08/2013 ver 2.12 by Makoto. Option 'amplimit' changed. Option 'lpf' help added.
% 10/08/2013 ver 2.11 by Makoto. Normalization computation debugged and refined. Plot titles for the case of combined conditions in the STUDY design fixed. %2.1f uniformly. Brushed up plot details.
% 10/04/2013 ver 2.10 by Makoto. Interaction (i.e. difference of difference) supported. Simplified difference plot title.
% 10/02/2013 ver 2.00 by Makoto. 6x faster after the 2nd run (31 vs 5 sec, 2200ICs). Stores precompute measure for maximum efficiency. Optimized the loop to generate topoerpconv_stack. Handles no-entry clusters correctly. Cleaned up warnings by auto code analyzer.
% 08/15/2013 ver 1.00 by Makoto. Version renumbered for official release.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 08/15/2013 ver 6.4 by Makoto. forcePositiveTopo added.
% 05/20/2013 ver 6.3 by Makoto. Bug STUDY.design(1,STUDY.currentdesign) fixed.
% 05/10/2013 ver 6.2 by Makoto. lpf supported. waitbar added. bsxfun used.
% 01/30/2013 ver 6.1 by Makoto. group_topo_now related bug fixed.
% 01/28/2013 ver 6.0 by Makoto. Combined groups supported. group_topo_now simplified.
% 10/10/2012 ver 5.3 by Makoto. 'sortedVariance' 'sortedVariancenorm' 'fillenvsumcolor' 'fillenvsum' added.
% 04/29/2012 ver 5.2 by Makoto. Bug fixed. STUDY.design(STUDY.currentdesign)
% 04/24/2012 ver 5.1 by Makoto. Revived 'amplimit' 'fillclust' options, and drawvertlineical doted lines for limitcontribtime range. Default changed into 'pvaf' from 'rv'. 
% 02/28/2012 ver 5.0 by Makoto. Values are now normalized by the number of subjects participating to the cluster. Group difference is better addressed when performing subtraction across groups. A group ratio is output too.
% 02/24/2012 ver 4.5 by Makoto. New color schema to improve color readability & shared line colors across multiple plots
% 02/22/2012 ver 4.4 by Makoto. Now plot title can accept numbers
% 01/24/2012 ver 4.3 by Makoto. Output added for the function.
% 01/17/2012 ver 4.2 by Makoto. Improved diff option; now mutually exclusive with normal plot. Diff title improved too. 
% 01/06/2012 ver 4.1 by Makoto. Time range options fixed. STUDY.design selection fixed. Wrong labels fixed. (Thanks to Agatha Lenartowicz!) 'sortedVariance' and 'topotype' implemented in the main function. Help edited.
% 12/30/2011 ver 4.0 by Makoto. Main function 100% redesigned. Full support for Study design.
% 12/26/2011 ver 3.2 by Makoto. Formats corrected.
% 12/25/2011 ver 3.1 by Makoto. Bug fixed about choosing wrong clusters under certain conditions. 
% 10/12/2011 ver 3.0 by Makoto. Completed the alpha version. 
% 10/11/2011 ver 2.2 by Makoto. Time range, group with no entry fixed
% 10/10/2011 ver 2.1 by Makoto. Two types of scalp topos can be presented.
% 10/08/2011 ver 2.0 by Makoto. Multiple groups supported (but one group one time using option). Result topos are now retrieved from the stored ones in STUDY 
% 08/01/2011 ver 1.0 by Makoto. Updated the original script to read data from STUDY. 180x faster than the original (3971 vs 22 sec, 700ICs).

% Copyright (C) 2016, Clement Lee, Makoto Miyakoshi, Hilit Serby, Arnold Delorme, Scott Makeig
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1.07  USA

function STUDY = std_envtopo(STUDY, ALLEEG, varargin)

% std_erpplot() with 'noplot' -> reads erp data and topography
% std_readtopo() reads topography
% std_readtopoclust

% if there is < 2 arguments, show help
if nargin < 2
    help std_envtopo;
    return
end

% if arguments are in odd number, wrong input
if mod(nargin,2) % if not an even number of arguments
    error('std_envtopo: Input argument list must be pairs of: ''keyx'', ''valx'' ');
end

% run finputcheck
arglist = finputcheck(varargin, ...
             {'clusterIdxToUse'      'integer'  -length(STUDY.cluster):length(STUDY.cluster)  [];...
              'conditions'           'integer'  1:length(STUDY.design(STUDY.currentdesign).variable(1,1).value) 1:length(STUDY.design(STUDY.currentdesign).variable(1,1).value);...
              'clust_exclude'        'integer'  1:length(STUDY.cluster)  [];...
              'onlyclus'             'string'   {'on', 'off'}  'on';...
              'statPvaf'             'real'     []  [];...
              'statPvafTail'         'integer'  []  [];...
              'baseline'             'real'     []  [];...
              'designVector'         'real'     []  [];...
              'timerange'            'real'     []  ALLEEG(1,1).times([1 end]);...
              'limitamp'             'real'     []  [];...
              'limitcontribtime'     'real'     []  [];...
              'drawvertline'         'real'     []  [];...
              'measure'              'string'   {'pvaf', 'ppaf', 'maxp', 'rltp', 'auc'}  'pvaf';...
              'normalizeClusterPvaf' 'string'   {'on', 'off'}  'off';...
              'fillenvsum'           'string'   {'selective','all','off'}  'selective';...
              'topotype'             'string'   {'inst', 'ave'}  'inst';...
              'fillenvsumcolor'      'real'     []  [0.875  0.875  0.875];...
              'fillclust'            'integer'  []  [];...
              'fillclustcolor'       'real'     []  [1 0 0];...
              'lpf'                  'real'     []  [];...
              'forcePositiveTopo'    'string'   {'on', 'off'}  'off';...
              'normalization'        'string'   {'normIC','normSubj','off'}  'normIC';...
              'downsampled'          'string'   {'on', 'off'}  'on';...
              });
          
if ischar(arglist)
    error(arglist);
end
clear varargin

% Sort clusterIdxToUse.
arglist.clusterIdxToUse = sort(arglist.clusterIdxToUse);

% Obtain erpTimes. Updated to be compatible with EEGLAB15. 08/10/2018 Makoto.
try % For EEGLAB version 14 or below.
    erpTimes = STUDY.cluster(1,2).erptimes;
catch % EEGLAB version 15 or later. Note this makes STUDY into classical mode, which should not be assign in to the base workspace.
    for clustIdx = 2:length(STUDY.cluster)
        [STUDY, datavals, erpTimes] = std_readdata(STUDY, ALLEEG, 'design', STUDY.currentdesign, 'clusters', clustIdx, 'datatype', 'erp', 'singletrials', 'off', 'timerange', [ALLEEG(1,1).xmin*1000 ALLEEG(1,1).xmax*1000]);
        STUDY.cluster(clustIdx).erpdata = datavals;
    end
    
    % If EEGLAB15 AND only one-way design, enter {''} to var2. This is a good place to also take care of this. 08/10/2018 Makoto.
    if length(STUDY.design(STUDY.currentdesign).variable) == 1
        STUDY.design(STUDY.currentdesign).variable(1,2).value = {''};
    end
end

% Calculate timepoints to plot and compute pvaf frames
arglist.totalPlotRange = find(erpTimes >= arglist.timerange(1)        & erpTimes <= arglist.timerange(2));
arglist.pvafFrames     = find(erpTimes >= arglist.limitcontribtime(1) & erpTimes <= arglist.limitcontribtime(2));

[~,arglist.relativeFramesWithinPlotRange] = intersect(arglist.totalPlotRange, arglist.pvafFrames);

% baseline period to timepoints (from sample)
if  ~isempty(arglist.baseline)
    arglist.baselinelimits = find(erpTimes >= arglist.baseline(1) & erpTimes <= arglist.baseline(2));
end

% Outlier cluster & all clusters included
try strcmp(STUDY.cluster(2).name(1:7), 'outlier');
    disp('Outlier cluster detected: Cls 3 and after will be used.')
    arglist.clust_all = 3:length(STUDY.cluster);% exclude parent (1) and outliers (2)
catch
    disp('No outlier cluster: Cls 2 and after will be used.')
    arglist.clust_all = 2:length(STUDY.cluster);% exclude parent (1) only
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Determine clusters to plot    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if arglist.clusterIdxToUse > 0
    arglist.clust_selected = arglist.clusterIdxToUse;
    arglist.clust_topentries = [];
elseif arglist.clusterIdxToUse < 0
    arglist.clust_selected = [];
    arglist.clust_topentries = abs(arglist.clusterIdxToUse);
end
if strcmp(arglist.onlyclus, 'on')
    arglist.clust_grandERP = setdiff(arglist.clust_all, arglist.clust_exclude);
else
    arglist.clust_grandERP = arglist.clust_all;
end
arglist.clustlabels = {STUDY.cluster.name};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Separate topos into groups   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(STUDY.group)>1 | length(STUDY.session)>1
    [STUDY, ALLEEG] = separateToposIntoGroups(STUDY, ALLEEG);
    % assignin('base', 'STUDY', STUDY); % Disabled because this does not respect EEGLAB15 changes.
    
    % Obtain group/session slot number and label. I believe this works with EEGLAB14 or below as well (untested!)
    variableLabels = {STUDY.design(STUDY.currentdesign).variable.label};
    isGroupVer1orVer2 = find(strcmp(lower(variableLabels), 'group') |  strcmp(lower(variableLabels), 'session'));
    if length(isGroupVer1orVer2)>1
        error('There is only one group or session supported.')
    end
    groupSessionLabel = lower(variableLabels{isGroupVer1orVer2});
else
    disp('No need to separate topos into groups.')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Convolve scalp topo with ERP   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numLoop = length(arglist.clust_grandERP)*...
    length(STUDY.design(1,STUDY.currentdesign).variable(1,1).value)*...
    length(STUDY.design(1,STUDY.currentdesign).variable(1,2).value);

convTopoErp = cell(1,length(arglist.clust_all));
if strcmp(arglist.downsampled,'on')
    originalConvTopoErpAllClusters = single(zeros(377, ALLEEG(1,1).pnts, length(arglist.clust_all), length(STUDY.design(1,STUDY.currentdesign).variable(1,1).value), length(STUDY.design(1,STUDY.currentdesign).variable(1,2).value)));
elseif strcmp(arglist.downsampled,'off')
    originalConvTopoErpAllClusters = single(zeros(3409, ALLEEG(1,1).pnts, length(arglist.clust_all), length(STUDY.design(1,STUDY.currentdesign).variable(1,1).value), length(STUDY.design(1,STUDY.currentdesign).variable(1,2).value)));
end

counter = 0;
waitbarHandle = waitbar(counter, 'Preparing the data. Please wait.');
for cls = 1:length(arglist.clust_all) % For all clusters for grandERP
    for columnwise = 1:size(STUDY.cluster(1,arglist.clust_all(cls)).erpdata, 1) % For the first variable
        for rowwise = 1:size(STUDY.cluster(1,arglist.clust_all(cls)).erpdata, 2) % For the second variable
            
            counter = counter+1;
            waitbar(counter/numLoop, waitbarHandle)
            
            if ~isempty(STUDY.cluster(1,arglist.clust_all(cls)).erpdata{columnwise, rowwise}) % Check if IC are present!
                currentClusterErps = single(STUDY.cluster(1,arglist.clust_all(cls)).erpdata{columnwise, rowwise}');
                
                % Low-pass filter the data if requested.
                if ~isempty(arglist.lpf)
                    currentClusterErps = myeegfilt(currentClusterErps, ALLEEG(1,1).srate,0,arglist.lpf);
                end
                
                % Obtain icawinv scalp topography.
                if  length(STUDY.group)>1 || length(STUDY.session)>1
                    if isGroupVer1orVer2 == 1
                        topo = STUDY.cluster(1,arglist.clust_all(cls)).topoall_group{columnwise, 1};
                    else
                        topo = STUDY.cluster(1,arglist.clust_all(cls)).topoall_group{1, rowwise};
                    end
                else
                    topo  = STUDY.cluster(1,arglist.clust_all(cls)).topoall;
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% convolve topo x erp %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                if strcmp(arglist.downsampled,'on')
                    tmpTopo0 = reshape(cell2mat(topo), [67 67 length(topo)]);
                    tmpTopo1 = tmpTopo0(1:3:67,1:3:67,:);
                    tmpTopo2 = tmpTopo1(~isnan(tmpTopo1));
                    tmpTopo3 = single(reshape(tmpTopo2, [length(tmpTopo2)/length(topo) length(topo)]));
                    convTopoErp{1, cls}{columnwise,rowwise} = tmpTopo3*currentClusterErps;  %min abs = zero
                    convTopoErpAllClusters(:,:,cls, columnwise, rowwise) = convTopoErp{1, cls}{columnwise,rowwise};
                    
                elseif strcmp(arglist.downsampled,'off')
                    tmpTopo1 = reshape(cell2mat(topo), [67 67 length(topo)]);
                    tmpTopo2 = tmpTopo1(~isnan(tmpTopo1));
                    tmpTopo3 = single(reshape(tmpTopo2, [length(tmpTopo2)/length(topo) length(topo)]));
                    convTopoErp{1, cls}{columnwise,rowwise} = tmpTopo3*currentClusterErps;  %min abs = zero
                    convTopoErpAllClusters(:,:,cls, columnwise, rowwise) = convTopoErp{1, cls}{columnwise,rowwise};
                end
                
            else
                disp([char(10) ' CAUTION: No IC entry in Cluster ' num2str(arglist.clust_all(1,1)+cls-1) ' Group ' num2str(rowwise)])
                convTopoErp{1, cls}{columnwise,rowwise} = single(zeros(size(convTopoErpAllClusters,1), ALLEEG(1,1).pnts)); % Dummy data
            end
            
            % topoerpconv_stack has 5 dimensions that are:
            % 1. topo*erp == 3409 or downsampled 377
            % 2. time points (which should be equal to ALLEEG(1,1).pnts
            % 3. cluster number
            % 4. variable 1 (usually within-subject  condition)
            % 5. variable 2 (usually between-subject condition)
        end
    end
end

close(waitbarHandle)
clear cls columnwise erp rowwise topo topoerpconv tmp* counter numLoop waitbarHandle
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Normalization (choice of by subj or by IC (default))%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This code is for EEGLAB14 or below.
if isfield(STUDY.cluster(1,2), 'setinds')
    
    setinds = {STUDY.cluster(arglist.clust_all).setinds};
    
    % If group is present, use this code.
    if length(STUDY.group)>1 || length(STUDY.session)>1
        numSubjectsByGroup = zeros(length(arglist.clust_all), length(STUDY.design(STUDY.currentdesign).variable(isGroupVer1orVer2).value));
        numIcsByGroup      = zeros(length(arglist.clust_all), length(STUDY.design(STUDY.currentdesign).variable(isGroupVer1orVer2).value));

        for clusterIdx = 1:length(arglist.clust_all)
            for groupIdx = 1:length(STUDY.design(STUDY.currentdesign).variable(isGroupVer1orVer2).value)
                if isGroupVer1orVer2 == 1
                    numIcsByGroup(clusterIdx,groupIdx)      = length(setinds{1,clusterIdx}{groupIdx,1});
                    numSubjectsByGroup(clusterIdx,groupIdx) = length(unique(setinds{1,clusterIdx}{groupIdx,1}));
                else
                    numIcsByGroup(clusterIdx,groupIdx)      = length(setinds{1,clusterIdx}{1,groupIdx});
                    numSubjectsByGroup(clusterIdx,groupIdx) = length(unique(setinds{1,clusterIdx}{1,groupIdx}));
                end
            end
        end
        
    % If group is not present, use this code.
    else
        numSubjectsByGroup = zeros(length(arglist.clust_all), 1);
        numIcsByGroup      = zeros(length(arglist.clust_all), 1);
        for clusterIdx = 1:length(arglist.clust_all)
            numIcsByGroup(clusterIdx)      = length(setinds{1,clusterIdx}{1,1});
            numSubjectsByGroup(clusterIdx) = length(unique(setinds{1,clusterIdx}{1,1}));
        end
        
    end
    
% This code is for EEGLAB15 or above.
else
    
    % Obtain group-mixed indices.
    setinds = {STUDY.cluster(arglist.clust_all).sets};
    
    % If between-subject condition is present.
    if length(STUDY.group)>1 || length(STUDY.session)>1
        
        % Identify unqiue group labels.
        if strcmp(groupSessionLabel, 'group')
            alleegGroupLabels = {ALLEEG.group};
        else
            alleegGroupLabels = {ALLEEG.session};
        end
        uniqueGroupLabels = unique(alleegGroupLabels);
        
        % Identify groupIdx vector.
        groupVector = zeros(length(alleegGroupLabels),1);
        for numUniqueLabelIdx = 1:length(uniqueGroupLabels)
            currentGroupSubjIdx = find(strcmp(alleegGroupLabels,uniqueGroupLabels{numUniqueLabelIdx}));
            for clsIdxIdx = 1:length(arglist.clust_all)
                currentClusterSets = setinds{clsIdxIdx};
                numIcsByGroup(clsIdxIdx,numUniqueLabelIdx) = length(find(ismember(currentClusterSets, currentGroupSubjIdx)));
                numSubjectsByGroup(clsIdxIdx,numUniqueLabelIdx) = length(find(ismember(unique(currentClusterSets), currentGroupSubjIdx)));
            end
        end
        
    % If between-subject condition is absent.
    else
        numIcsByGroup = cellfun('length', setinds)';
        uniqueSubjectsCell = cellfun(@(x) unique(x), setinds, 'uniformoutput', false);
        numSubjectsByGroup = cellfun('length', uniqueSubjectsCell)';
    end
end

% Perform normalization.
if strcmp(arglist.normalization, 'normSubj')
    normalizationFactors(1,1,1:size(numSubjectsByGroup,1),1,1:size(numSubjectsByGroup,2)) = numSubjectsByGroup;
elseif strcmp(arglist.normalization,'normIC')
    normalizationFactors(1,1,1:size(numSubjectsByGroup,1),1,1:size(numSubjectsByGroup,2)) = numIcsByGroup;
elseif strcmp(arglist.normalization, 'off')
    normalizationFactors(1,1,1:size(numSubjectsByGroup,1),1,1:size(numSubjectsByGroup,2)) = ones(size(numIcsByGroup));
end

% If group/session is present AND it is registered as var1.
if  length(STUDY.group)>1 || length(STUDY.session)>1
    if isGroupVer1orVer2 == 1
      normalizationFactors = permute(normalizationFactors, [1 2 3 5 4]);
    end
end

% Perform normalization.
convTopoErpAllClusters = bsxfun(@rdivide, convTopoErpAllClusters, normalizationFactors);
arglist.normalize_cluster = numSubjectsByGroup;



% This should work faster, but in reality it is 50-100x slower; why?
% invNsubjectGroup = (nsubject_group).^-1;
% normalizationFactors(1,1,1:size(invNsubjectGroup,1),1,1:size(invNsubjectGroup,2)) = invNsubjectGroup;
% convTopoErpAllClusters = bsxfun(@times, convTopoErpAllClusters, normalizationFactors);
% arglist.normalize_cluster = nsubject_group;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Report number of unique subjects for each cluster %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 1:size(arglist.normalize_cluster,1) % for the number of clusters
    str = sprintf('In %s, ', arglist.clustlabels{1,arglist.clust_all(n)});
    if length(STUDY.group)>1 % if multiple groups
        for m = 1:size(STUDY.group, 2)
            str = [str sprintf('%d %s ', arglist.normalize_cluster(n,m), num2str(STUDY.group{1,m}))]; %#ok<*AGROW>
            str = [str sprintf('(%d ICs) ', numIcsByGroup(n,m))]; %#ok<*AGROW>
        end
        disp(str)
    elseif length(STUDY.session)>1 % if multiple sessions
        for m = 1:size(STUDY.session, 2)
            str = [str sprintf('%d %s ', arglist.normalize_cluster(n,m), num2str(STUDY.session(1,m)))]; %#ok<*AGROW>
            str = [str sprintf('(%d ICs) ', numIcsByGroup(n,m))]; %#ok<*AGROW>
        end
        disp(str)
    else
        str = [str sprintf('%d %s ', arglist.normalize_cluster(n,1)) 'unique subjects'];
        str = [str sprintf('(%d ICs) ', numIcsByGroup(n,1))]; %#ok<*AGROW>            
        disp(str)
    end
end
clear column n row setinds nsubject* str tmp* whos* m userInput dataSizeIn* nIC_group

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Select cluseters to use %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,intersectIdx]       = intersect(arglist.clust_all, arglist.clust_grandERP);
convTopoErpAllClusters = convTopoErpAllClusters(:,:,intersectIdx,:,:);
arglist.clustlabels    = {STUDY.cluster(1,arglist.clust_grandERP).name}; % leave only the selected clusters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Apply new baseline (if specified) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
if isfield(arglist, 'baselinelimits')
    baselineValues = mean(convTopoErpAllClusters(:, arglist.baselinelimits,:,:,:), 2);
    convTopoErpAllClusters = bsxfun(@minus, convTopoErpAllClusters, baselineValues);
    fprintf('\nNew baseline applied.\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate outermost envelope for each conditions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
outermostEnv = sum(convTopoErpAllClusters, 3);
originalConvTopoErpAllClusters = convTopoErpAllClusters;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Select clusters (if specified) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if arglist.clusterIdxToUse > 0;
    [~,idx] = intersect(arglist.clust_grandERP, arglist.clusterIdxToUse);
    convTopoErpAllClusters = convTopoErpAllClusters(:,:,idx,:,:);
    tmp = {arglist.clustlabels(1, idx)};
    arglist.clustlabels = tmp{1,1}; clear tmp idx
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Process the conditions names into usable plot title. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var1Label = STUDY.design(STUDY.currentdesign).variable(1,1).value{1, arglist.designVector(1,1)};
var2Label = STUDY.design(STUDY.currentdesign).variable(1,2).value{1, arglist.designVector(1,2)};
var1Label = summarizeMultipleLabels(var1Label);
var2Label = summarizeMultipleLabels(var2Label);
if length(arglist.designVector)>2
    var3Label = STUDY.design(STUDY.currentdesign).variable(1,1).value{1, arglist.designVector(1,3)};
    var4Label = STUDY.design(STUDY.currentdesign).variable(1,2).value{1, arglist.designVector(1,4)};
    var3Label = summarizeMultipleLabels(var3Label);
    var4Label = summarizeMultipleLabels(var4Label);
    if length(arglist.designVector)>4
        var5Label = STUDY.design(STUDY.currentdesign).variable(1,1).value{1, arglist.designVector(1,5)};
        var6Label = STUDY.design(STUDY.currentdesign).variable(1,2).value{1, arglist.designVector(1,6)};
        var7Label = STUDY.design(STUDY.currentdesign).variable(1,1).value{1, arglist.designVector(1,7)};
        var8Label = STUDY.design(STUDY.currentdesign).variable(1,2).value{1, arglist.designVector(1,8)};
        var5Label = summarizeMultipleLabels(var5Label);
        var6Label = summarizeMultipleLabels(var6Label);
        var7Label = summarizeMultipleLabels(var7Label);
        var8Label = summarizeMultipleLabels(var8Label);
    end
end



%%%%%%%%%%%%%%%%%%%%
%%% Plot envtopo %%%
%%%%%%%%%%%%%%%%%%%%
figure
set(gcf,'Color', [0.93 0.96 1], 'NumberTitle', 'off', 'MenuBar', 'none', 'Name', 'std_envtopo()') % EEGLAB color
switch length(arglist.designVector)
    
    % Single condition.
    case 2
        arglist.title = [var1Label ' ' var2Label];
        envtopo_plot(STUDY, ALLEEG, squeeze(outermostEnv(:, arglist.totalPlotRange, 1, arglist.designVector(1,1), arglist.designVector(1,2))), squeeze(convTopoErpAllClusters(:,arglist.totalPlotRange,:,arglist.designVector(1,1),arglist.designVector(1,2))), arglist);

    % A-B.
    case 4
        arglist.title   = sprintf('{%s %s \n -%s %s}', var1Label, var2Label, var3Label, var4Label);
        outermostEnv    = outermostEnv(:,:,:,arglist.designVector(1,1),arglist.designVector(1,2)) - outermostEnv(:,:,:,arglist.designVector(1,3),arglist.designVector(1,4));
        diffEnv         = convTopoErpAllClusters(:,:,:,arglist.designVector(1,1),arglist.designVector(1,2))-convTopoErpAllClusters(:,:,:,arglist.designVector(1,3),arglist.designVector(1,4));
        [~, clusterIdx] = envtopo_plot(STUDY, ALLEEG, squeeze(outermostEnv(:,arglist.totalPlotRange,:,:,:)), squeeze(diffEnv(:,arglist.totalPlotRange,:,:,:)), arglist);

    % (A-B)-(C-D).
    case 8
        arglist.title   = sprintf('{%s %s \n -%s %s} \n-\n {%s %s \n -%s %s}', var1Label, var2Label, var3Label, var4Label, var5Label, var6Label, var7Label, var8Label);
        outermostEnv = outermostEnv(:,:,:,arglist.designVector(1,1),arglist.designVector(1,2)) - outermostEnv(:,:,:,arglist.designVector(1,3),arglist.designVector(1,4)) - (outermostEnv(:,:,:,arglist.designVector(1,5),arglist.designVector(1,6)) - outermostEnv(:,:,:,arglist.designVector(1,7),arglist.designVector(1,8)));
        diffEnv      = convTopoErpAllClusters(:,:,:,arglist.designVector(1,1),arglist.designVector(1,2)) - convTopoErpAllClusters(:,:,:,arglist.designVector(1,3),arglist.designVector(1,4)) - (convTopoErpAllClusters(:,:,:,arglist.designVector(1,5),arglist.designVector(1,6)) - convTopoErpAllClusters(:,:,:,arglist.designVector(1,7),arglist.designVector(1,8)));
        [~, clusterIdx] = envtopo_plot(STUDY, ALLEEG, squeeze(outermostEnv(:,arglist.totalPlotRange,:,:,:)), squeeze(diffEnv(:,arglist.totalPlotRange,:,:,:)), arglist);
end
        
  
                    
%%%%%%%%%%%%%%%%%%%%%
%%% Run statPvaf. %%%
%%%%%%%%%%%%%%%%%%%%%
if ~isempty(arglist.statPvaf)
    
    % 01/07/2020 Makoto. Updated.
    if arglist.clusterIdxToUse < 0
        error('To run PVAF stat, use only one cluster index.')
    end
    
    fprintf('Running statPvaf\n')

    % Obtain indices for Clusters
    if arglist.clusterIdxToUse > 0
        if strcmp(STUDY.cluster(2).name,'outlier 2')
            clusterIdx = arglist.clusterIdxToUse - 2; % adjust for parentcluster 1 and outlier 2
        else
            clusterIdx = arglist.clusterIdxToUse - 1; % adjust for parentcluster 1
        end
    end

    % Perform statistics
    statPvaf(STUDY, clusterIdx, arglist, originalConvTopoErpAllClusters, normalizationFactors);
end



% this function came with the original std_envtopo. I customized it only minimally. -Makoto
function [clusterPvafSorted, clusterIdx] = envtopo_plot(STUDY, ~, grandERP, projERP, varargin)

% data format conversion for compatibility
arglist = varargin{1,1};
clear varargin

for n = 1:size(projERP, 3);
    tmp{1,n} = projERP(:,:,n);
end
projERP = tmp;
clear n tmp

MAXTOPOS = 20;      % max topoplots to plot
drawvertlineWEIGHT = 2.0;  % lineweight of specified drawvertlineical lines
limitcontribtimeWEIGHT = 2.0; % lineweight of limonctrib drawvertlineical lines

currentFigureHandle = gcf; % remember the current figure (for Matlab 7.0.0 bug)

% add options (default values)
arglist.envmode     = 'avg';
arglist.voffsets    = [];
arglist.colorfile   = '';
arglist.colors      = '';
arglist.xlabel      = 'on';
arglist.ylabel      = 'on';
arglist.topoarg     = 0;

% empty cluster crashes the code, so 
for n = 2:length(STUDY.cluster)
    if ~isempty(STUDY.cluster(1,n).topoall)
        if strcmp(arglist.downsampled,'on')
            arglist.gridind = find(~isnan(STUDY.cluster(1,n).topoall{1,1}(1:3:67,1:3:67)));
        elseif strcmp(arglist.downsampled,'off')
            arglist.gridind = find(~isnan(STUDY.cluster(1,n).topoall{1,1}));
        end
        break
    end
end
clear n


frames = length(arglist.totalPlotRange);

if ~isempty(arglist.colors)
    arglist.colorfile = arglist.colors; % retain old usage 'colorfile' for 'colors' -sm 4/04
end

if ~isempty(arglist.drawvertline)
    arglist.drawvertline = arglist.drawvertline/1000; % condrawvertline from ms to s
end
%
%%%%%% Collect information about the gca, then delete it %%%%%%%%%%%%%
%
uraxes          = gca; % the original figure or subplot axes
currentPosition = get(uraxes,'Position');
axcolor         = get(uraxes,'Color');
delete(gca)

%
%%% Condrawvertline arglist.timerange and arglist.limitcontribtime to sec from ms %%%%
%
arglist.timerange = arglist.timerange/1000;   % the time range of the input data
arglist.limitcontribtime = arglist.limitcontribtime/1000; % the time range in which to select largest components


%
%%%%%%%%%%%% Collect time range information %%%%%%%%%%%%%%%%%%%%%%%%%%
%

xunitframes = 0;
xmin = arglist.timerange(1);
xmax = arglist.timerange(2);
pmin = xmin;
pmax = xmax;


%
%%%%%%%%%%%%%%% Find limits of the component selection window %%%%%%%%%
%

if arglist.limitcontribtime(1)<xmin
    arglist.limitcontribtime(1) = xmin;
end
if arglist.limitcontribtime(2)>xmax
    arglist.limitcontribtime(2) = xmax;
end
limframe1 = arglist.relativeFramesWithinPlotRange(1);
limframe2 = arglist.relativeFramesWithinPlotRange(end);

% Add dotted line to show pvaf frames
arglist.drawvertline(end+1) = arglist.limitcontribtime(1);
arglist.drawvertline(end+1) = arglist.limitcontribtime(2);



%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set 12 line colors.%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% 16 colors names officially supported by W3C specification for HTML
colors{1,1}  = [1 1 1];            % White
colors{2,1}  = [1 1 0];            % Yellow
colors{3,1}  = [1 0 1];            % Fuchsia
colors{4,1}  = [1 0 0];            % Red
colors{5,1}  = [0.75  0.75  0.75]; % Silver
colors{6,1}  = [0.5 0.5 0.5];      % Gray
colors{7,1}  = [0.5 0.5 0];        % Olive
colors{8,1}  = [0.5 0 0.5];        % Purple
colors{9,1}  = [0.5 0 0];          % Maroon
colors{10,1} = [0 1 1];            % Aqua
colors{11,1} = [0 1 0];            % Lime
colors{12,1} = [0 0.5 0.5];        % Teal
colors{13,1} = [0 0.5 0];          % Green
colors{14,1} = [0 0 1];            % Blue
colors{15,1} = [0 0 0.5];          % Navy
colors{16,1} = [0 0 0];            % Black

% Silver is twice brighter because used for background
colors{5,1} = [0.875 0.875 0.875];

% Choosing and sorting 12 colors for line plot, namely Red, Blue, Green, Fuchsia, Lime, Aqua, Maroon, Olive, Purple, Teal, Navy, and Gray
linecolors = colors([4 13 14 3 11 10 9 7 8 12 15 6]);



%
%%%%%%%%%%%%%%%% Check other input variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
arglist.clusterIdxToUseOriginal = arglist.clusterIdxToUse;
[tmp,pframes] = size(projERP{1});
if frames ~= pframes
    error('Size of trial in projERP and grandenvERP do not agree');
end

if isempty(arglist.voffsets) || ( size(arglist.voffsets) == 1 && arglist.voffsets(1) == 0 )
    arglist.voffsets = zeros(1,MAXTOPOS);
end

if isempty(arglist.clusterIdxToUse) || arglist.clusterIdxToUse(1) == 0
   arglist.clusterIdxToUse = 1:length(projERP); % by default the number of projected ERP input
end
if min(arglist.clusterIdxToUse) < 0
    if length(arglist.clusterIdxToUse) > 1
        error('Negative clustnums must be a single integer.');
    end
    if -arglist.clusterIdxToUse > MAXTOPOS
        fprintf('Can only plot a maximum of %d components.\n',MAXTOPOS);
        return
    else
        MAXTOPOS = -arglist.clusterIdxToUse;
        arglist.clusterIdxToUse = 1:length(projERP); % by default the number of projected ERP input
    end
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Process components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
numClusters = length(arglist.clusterIdxToUse);

%
%%%%%%%%%%%%%%% Compute plotframes and envdata %%%%%%%%%%%%%%%%%%%%%
%
numTopos = length(arglist.clusterIdxToUse);
if numTopos > MAXTOPOS
    numTopos = MAXTOPOS; % limit the number of topoplots to display
end


%
% first, plot the data envelope
%
envdata = zeros(2,frames*(numClusters+1));
envdata(:,1:frames) = envelope(grandERP, arglist.envmode);

fprintf('Data epoch is from %.0f ms to %.0f ms.\n',1000*xmin,1000*xmax);
fprintf('Plotting data from %.0f ms to %.0f ms.\n',1000*xmin,1000*xmax);
fprintf('Comparing maximum projections for components:  ');
if numClusters>32
    fprintf('\n');
end
compvars = zeros(1,numClusters);

%
% Compute frames to plot
%
sampint  = (xmax-xmin)/(frames-1);     % sampling interval in sec
times    = xmin:sampint:xmax;          % make vector of data time values

[~, minf] = min(abs(times-pmin));
[~, maxf] = min(abs(times-pmax));
if limframe1 < minf
    limframe1 = minf;
end
if limframe2 > maxf
    limframe2 = maxf;
end

%
%%%%%%%%%%%%%% find max variances and their frame indices %%%%%%%%%%%
%
pvafPeakFrames = zeros(numClusters,1);
for clusterIdx = 1:numClusters
    if ~rem(clusterIdx,5)
        fprintf('%d ... ',arglist.clusterIdxToUse(clusterIdx)); % clusterIdx is index into clustnums
    end
    if ~rem(clusterIdx,100)
        fprintf('\n');
    end

    envdata(:,clusterIdx*frames+1:(clusterIdx+1)*frames) = envelope(projERP{clusterIdx}, arglist.envmode);

    [maxval,maxi] = max(sum(projERP{clusterIdx}(:,limframe1:limframe2).*projERP{clusterIdx}(:,limframe1:limframe2)));
    % find point of max variance for comp clusterIdx
    compvars(clusterIdx)   = maxval;
    maxi = maxi+limframe1-1;
    pvafPeakFrames(clusterIdx) = maxi;
    maxproj(:,clusterIdx)  = projERP{clusterIdx}(:,maxi); % Note: maxproj contains only arglist.plotchans -sm 11/04

end % component clusterIdx
fprintf('\n');

%
%%%%%%%%%%%%%%% Compute component selection criterion %%%%%%%%%%%%%%%%%%%%%%%%%%
%

if ~xunitframes
    fprintf('  in the interval %3.0f ms to %3.0f ms.\n',1000*times(limframe1),1000*times(limframe2));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate all-cluster PVAF. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clusterPvaf  = zeros(length(projERP),1);
allClusterProjection = grandERP(:,limframe1:limframe2);
for clustIdx = 1:length(projERP)
    
    remainingClusterProjection = grandERP(:,limframe1:limframe2)-projERP{clustIdx}(:,limframe1:limframe2);
    
    if strcmpi(arglist.measure,'pvaf')     % Percent Variance Accounted For (PVAF).
        clusterPvaf(clustIdx) = 100-100*mean(var(remainingClusterProjection))/mean(var(allClusterProjection));
        fprintf('%s variance %2.4f, percent variance accounted for (PVAF): %2.1f%%\n',arglist.clustlabels{clustIdx}, mean(var(projERP{clustIdx}(:,limframe1:limframe2))), clusterPvaf(clustIdx));
    
    elseif strcmpi(arglist.measure,'ppaf') % Percent Power Accounted For (PPAF).
        clusterPvaf(clustIdx) = 100-100*mean(mean(remainingClusterProjection.^2))/mean(mean(allClusterProjection.^2));
        fprintf('%s power %2.4f, percent power accounted for (PPAF): %2.1f%%\n',      arglist.clustlabels{clustIdx}, mean(mean(projERP{clustIdx}(:,limframe1:limframe2).^2)), clusterPvaf(clustIdx));
        
    elseif     strcmpi(arglist.measure,'maxp')  % Maximum Power of backproj
        clusterPvaf(clustIdx) = max(mean(projERP{clustIdx}(:,limframe1:limframe2).*projERP{clustIdx}(:,limframe1:limframe2)));
        fprintf('%s maximum mean power of back-projection: %g\n',arglist.clustlabels{clustIdx},clusterPvaf(clustIdx));
   
    elseif strcmpi(arglist.measure,'rltp')    % Relative Power
        powdat = mean(mean(grandERP(:,limframe1:limframe2).^2));
        clusterPvaf(clustIdx) = 100*mean(mean((projERP{clustIdx}(:,limframe1:limframe2)).^2))/powdat;
        fprintf('%s relative power of back-projection: %2.1f%%\n',arglist.clustlabels{clustIdx},clusterPvaf(clustIdx));
        
    elseif strcmpi(arglist.measure,'auc')    % Enveloped Area ratio
        outermostEnv      = envelope(grandERP(:,limframe1:limframe2), arglist.envmode);
        outermostEnv_auc = sum(abs(outermostEnv(1,:)))+sum(abs(outermostEnv(2,:)));
        tmp_projERP       = envelope(projERP{clustIdx}(:,limframe1:limframe2), arglist.envmode);
        tmp_projERP_auc  = sum(abs(tmp_projERP(1,:)))+sum(abs(tmp_projERP(2,:)));
        clusterPvaf(clustIdx)        = 100*tmp_projERP_auc/outermostEnv_auc;
        fprintf('%s enveloped auc ratio : %2.1f%%\n',arglist.clustlabels{clustIdx},clusterPvaf(clustIdx));
    end
end

%%%%%% Normalize the clusterPvaf to 100 if requested %%%%%%
if strcmp(arglist.normalizeClusterPvaf, 'on')
    clusterPvaf = clusterPvaf*1/sum(clusterPvaf)*100;
    fprintf('\n')
    fprintf('normalizeClusterPvaf is on: clusterPvaf is normalized to 100%% in the plot \n');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Sort clusters by PVAF. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[clusterPvafSorted, pvafSortingIdx] = sort(clusterPvaf, 'descend');
clusterIdxToUseSortedByPvaf         = arglist.clusterIdxToUse(pvafSortingIdx); % Actual cluster Index -1 (if no outlier cluster in Cls2) or -2 (if outlier cluster in Cls2)
pvafPeakLatencySortedByPvaf         = pvafPeakFrames(         pvafSortingIdx);              
maxprojSrotedByPvaf                 = maxproj(:,              pvafSortingIdx);         
clustLabelsSortedByPvaf             = arglist.clustlabels(    pvafSortingIdx);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Select clusters from the top. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clusterIdxToUseSortedAndSelected = clusterIdxToUseSortedByPvaf(1:numTopos);
clusterPvafSortedAndSelected     = clusterPvafSorted(          1:numTopos);
pvafPeakFramesSortedAndSelected  = pvafPeakLatencySortedByPvaf(1:numTopos);
maxprojSortedAndSelected         = maxprojSrotedByPvaf(:,      1:numTopos);
clustLabelsSortedAndSelected     = clustLabelsSortedByPvaf(    1:numTopos);% actual component numbers (output var)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Sort the select clusters by peak latency. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[selectedPeakLatencySortedByLatency, pvafPeakLatencySortIdx] = sort(pvafPeakFramesSortedAndSelected); % sort plotframes on their temporal order (min:max)
selectedClusterIdxSortexByLatency    = clusterIdxToUseSortedAndSelected(pvafPeakLatencySortIdx);
selectedClusterPvafSortedByLatency   = clusterPvafSortedAndSelected(    pvafPeakLatencySortIdx);
selectedMaxprojSortedByLatency       = maxprojSortedAndSelected(      :,pvafPeakLatencySortIdx);
selectedClusterLabelsSortedByLatency = clustLabelsSortedAndSelected(    pvafPeakLatencySortIdx);



% Set up graph axes.
vlen = length(arglist.voffsets); % extend voffsets if necessary
while vlen< numTopos
    arglist.voffsets = [arglist.voffsets arglist.voffsets(vlen)]; % repeat last offset given
    vlen=vlen+1;
end

head_sep = 1.2;
topowidth = currentPosition(3)/(numTopos+(numTopos-1)/5); % width of each topoplot
if topowidth > 0.20    % adjust for maximum height
    topowidth = 0.2;
end

if rem(numTopos,2) == 1  % odd number of topos
    topoleft = currentPosition(3)/2 - (floor(numTopos/2)*head_sep + 0.5)*topowidth;
else % even number of topos
    topoleft = currentPosition(3)/2 - (floor(numTopos/2)*head_sep)*topowidth;
end

%
%%%%%%%%%%%%%%%%%%%% Print times and frames of comp maxes %%%%%%%%%%%%%%
%
peakLatencyInS = times(selectedPeakLatencySortedByLatency);
if ~xunitframes
    fprintf('    with max var at times (ms): ');
    for t=1:numTopos
        fprintf('%4.0f  ',1000*peakLatencyInS(t));
    end
    fprintf('\n');
end

fprintf('                  epoch frames: ');
for clusterIdx = 1:numTopos
    fprintf('%4d  ',limframe1-1+selectedPeakLatencySortedByLatency(clusterIdx));
end
fprintf('\n');
fprintf('IC cluster %s in interval:  ', upper(arglist.measure));
fprintf('%2.1f ', selectedClusterPvafSortedByLatency);
fprintf('\n');

%
%%%%%%%%%%%%%%%%%%%%% Plot the data envelopes %%%%%%%%%%%%%%%%%%%%%%%%%
%
BACKCOLOR = [0.7 0.7 0.7];
newaxes=axes('position',currentPosition);
axis off
set(newaxes,'FontSize',16,'FontWeight','Bold','Visible','off');
set(newaxes,'Color',BACKCOLOR); % set the background color
delete(newaxes) %XXX

% site the plot at bottom of the current axes
axe = axes('Position',[currentPosition(1) currentPosition(2) currentPosition(3) 0.6*currentPosition(4)],...
    'FontSize',16,'FontWeight','Bold');

set(axe,'GridLineStyle',':')
set(axe,'Xgrid','off')
set(axe,'Ygrid','on')
axes(axe)
set(axe,'Color',axcolor);

% Plot the envelope of the summed selected components.
sumOfAllSelectedClusterProjections = zeros(size(projERP{1}));
for n = 1:numTopos
    if sum(arglist.clusterIdxToUseOriginal) < 0
        sumOfAllSelectedClusterProjections = sumOfAllSelectedClusterProjections + projERP{selectedClusterIdxSortexByLatency(n)}; % add up all cluster projections
    else
        sumOfAllSelectedClusterProjections = sumOfAllSelectedClusterProjections + projERP{pvafSortingIdx(n)}; % add up all cluster projections
    end
end
remainingClusterProjection = allClusterProjection-sumOfAllSelectedClusterProjections(:,limframe1:limframe2);

% Calculate PVAF of all the selected clusters.
if     strcmpi(arglist.measure,'pvaf')    % Percent Variance
    allSelectedClusterPvaf = 100-100*mean(var(remainingClusterProjection))/mean(var(allClusterProjection));

elseif strcmpi(arglist.measure,'ppaf') % Percent Power
    allSelectedClusterPvaf = 100-100*mean(mean(remainingClusterProjection.^2))/mean(mean(allClusterProjection.^2));
            
elseif strcmpi(arglist.measure,'maxp')  % Maximum Power of backproj
    allSelectedClusterPvaf = max(mean(remainingClusterProjection.*remainingClusterProjection));
    
elseif strcmpi(arglist.measure,'rltp')    % Relative Power
    powdat = mean(mean(grandERP(:,limframe1:limframe2).^2));
    allSelectedClusterPvaf = 100*mean(mean((remainingClusterProjection).^2))/powdat;

elseif strcmpi(arglist.measure,'auc')    % Enveloped Area ratio
    outermostEnv     = envelope(grandERP(:,limframe1:limframe2), arglist.envmode);
    outermostEnv_auc = sum(abs(outermostEnv(1,:)))+sum(abs(outermostEnv(2,:)));
    tmp_projERP      = envelope(remainingClusterProjection, arglist.envmode);
    tmp_projERP_auc  = sum(abs(tmp_projERP(1,:)))+sum(abs(tmp_projERP(2,:)));
    allSelectedClusterPvaf = 100*tmp_projERP_auc/outermostEnv_auc;
end

fprintf('Summed cluster %s in interval [%s %s] ms: %3.1f%', upper(arglist.measure), int2str(1000*times(limframe1)), int2str(1000*times(limframe2)), allSelectedClusterPvaf);
fprintf('\n')

%
% Plot the summed projection filled
%
outermostEnvelope = envelope(sumOfAllSelectedClusterProjections, arglist.envmode);
if     strcmpi(arglist.fillenvsum,'selective')
    filltimes = times(limframe1:limframe2);
    mins = matsel(outermostEnvelope,frames,0,2,0);
    p=fill([filltimes filltimes(end:-1:1)],...
         [matsel(outermostEnvelope,frames,limframe1:limframe2,1,0) mins(limframe2:-1:limframe1)], arglist.fillenvsumcolor);
    set(p,'EdgeColor', arglist.fillenvsumcolor);
    hold on
    % redraw outlines
    p=plot(times,matsel(envdata,frames,0,1,1),'Color', colors{16,1});% plot the max
    set(p,'LineWidth',2);
    p=plot(times,matsel(envdata,frames,0,2,1),'Color', colors{16,1});% plot the min
    set(p,'LineWidth',2);
    
elseif strcmpi(arglist.fillenvsum,'all')
    mins = matsel(outermostEnvelope,frames,0,2,0);
    p=fill([times times(frames:-1:1)],...
        [matsel(outermostEnvelope,frames,0,1,0) mins(frames:-1:1)], arglist.fillenvsumcolor);
    set(p,'EdgeColor', arglist.fillenvsumcolor);
    hold on
    % redraw outlines
    p=plot(times,matsel(envdata,frames,0,1,1),'Color', colors{16,1});% plot the max
    set(p,'LineWidth',2);
    p=plot(times,matsel(envdata,frames,0,2,1),'Color', colors{16,1});% plot the min
    set(p,'LineWidth',2);

else % if no 'fill'
    tmp = matsel(outermostEnvelope,frames,0,2,0);
    p=plot(times,tmp);% plot the min
    hold on
    set(p,'Color', colors{5,1});
    set(p,'linewidth',2);
    p=plot(times,matsel(outermostEnvelope,frames,0,1,0));% plot the max
    set(p,'linewidth',2);
    set(p,'Color', colors{5,1});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot the envelopes. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assign unique combinations of colors and line styles to a cluster. 08/13/2018 Makoto. 
lineStyleList = {'-' '--' ':'  '-.'};

% Plot the outermost envelope.
hold on
upperEnv = matsel(envdata,frames,0,1,1);
currentPlotHandle = plot(times,upperEnv);
set(currentPlotHandle, 'Color', colors{16,1},'LineWidth', 2, 'LineStyle', '-');
lowerEnv = matsel(envdata,frames,0,2,1);
currentPlotHandle = plot(times,lowerEnv);
set(currentPlotHandle, 'Color', colors{16,1},'LineWidth', 2, 'LineStyle', '-');

% Plot all inner envelopes.
if sum(arglist.clusterIdxToUseOriginal) < 0 % For the case of taking top n clusters.
    clusterSortingIdxByPeakLatency = arglist.clust_all(selectedClusterIdxSortexByLatency);
else % For the case of taking specified n clusters
    [~,sortingIdx]                     = sort(selectedClusterIdxSortexByLatency);
    [~,clusterSortingIdxByPeakLatency] = sort(sortingIdx);
end

for envelopeIdx = 1:numTopos

    % Set line color.
    if sum(arglist.clusterIdxToUseOriginal) < 0 % For the case of taking top n clusters.
        rankInTheAllClusters = find(arglist.clust_all==clusterSortingIdxByPeakLatency(envelopeIdx));
        currentMatselIdx     = rankInTheAllClusters+1;
    else
        rankInTheAllClusters = find(arglist.clust_all==selectedClusterIdxSortexByLatency(envelopeIdx));
        currentMatselIdx     = clusterSortingIdxByPeakLatency(envelopeIdx)+1;
    end
    
    lineColorIdx = mod(rankInTheAllClusters, size(linecolors,1)); % envelopeIdx-1 is an index for the original IC index.
    if lineColorIdx == 0
        lineColorIdx = size(linecolors,1);
    end
    
    % Set line style.
    lineStyleIdx = 1+floor(rankInTheAllClusters/size(linecolors,1));
    
    % Plot upper envelope.
    currentEnvelope = matsel(envdata,frames,0,1,currentMatselIdx); % +1 is because the first one is the outermost envelope.

    % Plot the envelope.
    currentPlotHandle = plot(times,currentEnvelope);
    set(currentPlotHandle, 'Color', linecolors{lineColorIdx}, 'LineStyle', lineStyleList{lineStyleIdx}, 'LineWidth', 1);

    % Plot lower envelope.
    currentEnvelope = matsel(envdata,frames,0,2,currentMatselIdx);
    currentPlotHandle = plot(times,currentEnvelope);
    set(currentPlotHandle, 'Color', linecolors{lineColorIdx}, 'LineStyle', lineStyleList{lineStyleIdx}, 'LineWidth', 1);
end
set(gca,'FontSize',14)



%%%%%%%%%%%%%%%%%%%%%%%%
%%% axis min and max %%%
%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(arglist.limitamp)
    ymin = min([grandERP(:); outermostEnvelope(:)]);
    ymax = max([grandERP(:); outermostEnvelope(:)]);
    yspace = (ymax-ymin)*0.01;
    ymin = ymin-yspace;
    ymax = ymax+yspace;
else
    ymin = arglist.limitamp(1);
    ymax = arglist.limitamp(2);
    ylim([ymin ymax])
end
xlim([pmin pmax])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% fill specified component %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(arglist.fillclust)
    if      sum(arglist.clusterIdxToUseOriginal) < 0
        [~,fillIdx] = ismember(arglist.fillclust, arglist.clust_grandERP(selectedClusterIdxSortexByLatency));
    elseif  sum(arglist.clusterIdxToUseOriginal) > 0
        [~,fillIdx] = ismember(arglist.fillclust, arglist.clust_selected);
    else
           fillIdx  = 0;
    end
        
    if any(fillIdx)
        fprintf('filling the envelope of component %d\n',arglist.fillclust);
        
        % Obtain matsel input.
        if sum(arglist.clusterIdxToUseOriginal) < 0
            currentFillclustIdx = find(arglist.clust_all==arglist.fillclust);
        else
            currentFillclustIdx = find(arglist.clusterIdxToUse==arglist.fillclust);
        end
        currentMatselIdx    = currentFillclustIdx+1;
        
        % Obtain the envelope.
        upperEnv = matsel(envdata,frames,0,1,currentMatselIdx);
        lowerEnv = matsel(envdata,frames,0,2,currentMatselIdx);
        
        % Plot fill.
        fill([times times(frames:-1:1)], [upperEnv lowerEnv(frames:-1:1)], arglist.fillclustcolor);
        
        % Overplot the data envlope again so it is not covered by the filled component
        plot(times,matsel(envdata,frames,0,1,1), 'k', 'LineWidth',2);% plot the max
        plot(times,matsel(envdata,frames,0,2,1), 'k', 'LineWidth',2);% plot the min
    else
        fprintf('cluster %d is not on the list for plotting\n',arglist.fillclust);
    end
else
    fillIdx = 0;
end
clear clustlist clustn



%%%%%%%%%%%%%%%%%%%%%%%%
%%% Add drawvertline %%%
%%%%%%%%%%%%%%%%%%%%%%%%  
% draw a vertlcal line at time zero
if arglist.timerange(1) <= 0 && arglist.timerange(2) >= 0
    vl=plot([0 0], [-1e10 1e10],'k');
    set(vl,'linewidth',1,'linestyle','-');
end

% if specified by option, draw a vertlcal line at specified latency
if ~isempty(arglist.drawvertline)
    for v=1:length(arglist.drawvertline)
        vl=plot([arglist.drawvertline(v) arglist.drawvertline(v)], [-1e10 1e10],'k');
        if any(arglist.limitcontribtime ~= 0) && v>= length(arglist.drawvertline)-1;
            set(vl,'linewidth',limitcontribtimeWEIGHT);
            set(vl,'linestyle',':');
        else
            set(vl,'linewidth',drawvertlineWEIGHT);
            set(vl,'linestyle','--');
        end
    end
end
ylim([ymin ymax])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Show clusterPvafSorted and fillclust %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = text(double(xmin+0.02*(xmax-xmin)), double(ymin+0.12*(ymax-ymin)), ...
    ['Summed ' upper(arglist.measure) ' ' num2str(allSelectedClusterPvaf,'%2.1f')]);
set(t,'fontsize',15)

if ~isempty(arglist.fillclust) && any(fillIdx)
    t = text(double(xmin+0.02*(xmax-xmin)), ...
        double(ymin+0.03*(ymax-ymin)), ...
        ['Cluster ' num2str(arglist.fillclust) ' filled.']);
    set(t,'fontsize',15)
end



%
%%%%%%%%%%%%%%%%%%%%%% Label axes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
set(axe,'Color',axcolor);

if strcmpi(arglist.xlabel, 'on')
    l= xlabel('Latency (s)');
    set(l,'FontSize',14);
end

if strcmpi(arglist.ylabel, 'on')
    l = ylabel('RMS uV');
    set(l,'FontSize',14);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Draw oblique lines. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
axall = axes('Position', currentPosition, 'Visible','Off','Fontsize',16); % whole-figure invisible axes
axes(axall)
set(axall,'Color',axcolor);
axis([0 1 0 1])
width  = xmax-xmin;
pwidth = pmax-pmin;
height = ymax-ymin;
hold(axall)
for topoIdx = 1:numTopos

    % Set line color.
    if sum(arglist.clusterIdxToUseOriginal) < 0 % For the case of taking top n clusters.
        rankInTheAllClusters = find(arglist.clust_all==clusterSortingIdxByPeakLatency(topoIdx));
        currentMatselIdx     = rankInTheAllClusters+1;
    else
        rankInTheAllClusters = find(arglist.clust_all==selectedClusterIdxSortexByLatency(topoIdx)); % This assigns fixed colors for cluster indices.
        currentClusterIdx    = clusterSortingIdxByPeakLatency(topoIdx);
        currentMatselIdx     = currentClusterIdx+1;
    end
    
    lineColorIdx = mod(rankInTheAllClusters, size(linecolors,1)); % envelopeIdx-1 is an index for the original IC index.
    if lineColorIdx == 0
        lineColorIdx = size(linecolors,1);
    end
    
    % Set line style.
    lineStyleIdx = 1+floor(rankInTheAllClusters/size(linecolors,1));
    
    % Obtain max envelope value.
    maxenv = matsel(envdata,frames,selectedPeakLatencySortedByLatency(topoIdx),1,currentMatselIdx);

    % Plot the oblique line and obtain the handle.
    data_y = 0.6*(arglist.voffsets(topoIdx)+maxenv-ymin)/height;
    if (data_y > currentPosition(2)+0.6*currentPosition(4))
        data_y = currentPosition(2)+0.6*currentPosition(4);
    end
    currentObliqueLineHandle = plot([(peakLatencyInS(topoIdx)-pmin)/pwidth topoleft + 1/currentPosition(3)*(topoIdx-1)*1.2*topowidth + (topowidth*0.6)], [data_y 0.653]); % 0.68 is bottom of topo maps
    set(currentObliqueLineHandle, 'Color', linecolors{lineColorIdx}, 'LineStyle', lineStyleList{lineStyleIdx}, 'LineWidth', 1);
end



%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot topoplots. %%%
%%%%%%%%%%%%%%%%%%%%%%%
for topoIdx = 1:numTopos % Left to right.
    currentTopoplotHandle = axes('Units','Normalized','Position',...
        [currentPosition(3)*topoleft+currentPosition(1)+(topoIdx-1)*head_sep*topowidth currentPosition(2)+0.66*currentPosition(4) ...
        topowidth topowidth*head_sep]);
    axes(currentTopoplotHandle);                         

    if strcmp(arglist.downsampled, 'on')
        if arglist.gridind ~= 0
            currentTopo = zeros(23,23);
            currentTopo(:)=nan ;
            currentTopo(arglist.gridind) = selectedMaxprojSortedByLatency(:,topoIdx);
        end
    elseif strcmp(arglist.downsampled, 'off')
        if arglist.gridind ~= 0
            currentTopo = zeros(67,67);
            currentTopo(:)=nan ;
            currentTopo(arglist.gridind) = selectedMaxprojSortedByLatency(:,topoIdx);
        end
    end

    % If arglist.forcePositiveTopo is 'on', flip the negative topos. 08/15/2013 Makoto
    if strcmp(arglist.forcePositiveTopo, 'on')
        % check topo polarity
        [~,oneDaddress] = max(abs(currentTopo(:)));
        [row,column] = ind2sub(size(currentTopo),oneDaddress);
        % if absolute maximum is negative, multiply -1 to reverse negative topo to positive 
        if currentTopo(row,column)<0
            currentTopo = currentTopo*-1;
        end
    end

    figure(currentFigureHandle);
    if strcmp(arglist.topotype, 'inst')
        toporeplot(currentTopo, 'style', 'both', 'plotrad',0.5,'intrad',0.5, 'verbose', 'off');
    else % which is arglist.topytype == 'ave'
        toporeplot(STUDY.cluster(1, arglist.clust_grandERP(selectedClusterIdxSortexByLatency(topoIdx))).topo, 'style', 'both', 'plotrad',0.5,'intrad',0.5, 'verbose', 'off');
    end

    % Label components.
    currentComplabel = selectedClusterLabelsSortedByLatency(topoIdx);
    text(0.00,  0.80, currentComplabel,'FontSize',14, 'FontWeight','Bold','HorizontalAlignment','Center');
    text(0.00, -0.605,[upper(arglist.measure) ' ' sprintf('%2.1f', selectedClusterPvafSortedByLatency(topoIdx))],'FontSize',13, 'HorizontalAlignment','Center');
end

% Use the colormap 'jet'.
colormap('jet')

% Plot a colorbar.
currentTopoplotHandle = axes('Position',[currentPosition(1)+currentPosition(3)*1.01  currentPosition(2)+0.675*currentPosition(4) currentPosition(3)*.02 currentPosition(4)*0.09]);
h=cbar(currentTopoplotHandle);   
set(h,'YAxisLocation', 'left', 'YTick', [0 1], 'YTickLabel', {'-' '+'}, 'FontSize', 16)
axes(axall)
set(axall,'Color',axcolor);



%%%%%%%%%%%%%%%%%%%%%%%
%%% plot main title %%%
%%%%%%%%%%%%%%%%%%%%%%%
%tmp = text(0.50,0.92, arglist.title,'FontSize',18,'HorizontalAlignment','Center','FontWeight','Bold');
tmp = text(0.50,1.01, arglist.title,'FontSize',14,'HorizontalAlignment','Center','FontWeight','Bold'); % 09/18/2018 Makoto.
set(tmp, 'interpreter', 'none');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% bring component lines to top %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(axall)
set(axall,'layer','top');



%
%%%%%%%%%%%%%%%%%%%%%%%%% turn on axcopy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
axcopy(gcf, 'if ~isempty(get(gca, ''''userdata'''')), eval(get(gca, ''''userdata'''')); end;');

return %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function envdata = envelope(data, envmode)  % also in release as env()
if nargin < 2
    envmode = 'avg';
end;
if strcmpi(envmode, 'rms');
    warning off; %#ok<*WNOFF>
    negflag = (data < 0);
    dataneg = negflaarglist.* data;
    dataneg = -sqrt(sum(datanearglist.*dataneg,1) ./ sum(negflag,1));
    posflag = (data > 0);
    datapos = posflaarglist.* data;
    datapos = sqrt(sum(datapos.*datapos,1) ./ sum(posflag,1));
    envdata = [datapos;dataneg];
    warning on; %#ok<*WNON>
else
    if size(data,1)>1
        maxdata = max(data); % max at each time point
        mindata = min(data); % min at each time point
        envdata = [maxdata;mindata];
    else
        maxdata = max([data;data]); % max at each time point
        mindata = min([data;data]); % min at each time point
        envdata = [maxdata;mindata];
    end
end;

return %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [STUDY, ALLEEG] = separateToposIntoGroups(STUDY, ALLEEG)

% Obtain group/session slot index and label.
variableLabels  = {STUDY.design(STUDY.currentdesign).variable.label};
isGroupVer1orVer2 = find(strcmp(lower(variableLabels), 'group') |  strcmp(lower(variableLabels), 'session'));
if length(isGroupVer1orVer2)>1
    error('There is only one group or session supported.')
end
groupSessionLabel = lower(variableLabels{isGroupVer1orVer2});

if strcmp(groupSessionLabel, 'group')
    alleegGroupLabels = {ALLEEG.group};
else
    alleegGroupLabels = {ALLEEG.session};
end

% Check if any item is numeric. 09/17/2018 Makoto.
for loopIdx = 1:length(alleegGroupLabels)
    if isnumeric(alleegGroupLabels{loopIdx})
        alleegGroupLabels{loopIdx} = num2str(alleegGroupLabels{loopIdx});
    end
end

uniqueGroupLabels = unique(alleegGroupLabels);
    
% This is the case for EEGLAB14 or below.
if isfield(STUDY.design(STUDY.currentdesign), 'cell') 
    numGroups = size(STUDY.design(STUDY.currentdesign).variable(isGroupVer1orVer2).value, 2);
    for clsIdx = 2:length(STUDY.cluster)
        
        % Delete existing topoall_group.
        if isfield(STUDY.cluster(1,clsIdx), 'topoall_group')
            STUDY.cluster(1,clsIdx).topoall_group = [];
        end
        
        for groupIdx = 1:numGroups
            if isGroupVer1orVer2 == 1
                for icIdx = 1:length(STUDY.cluster(1,clsIdx).setinds{groupIdx,1})
                    tmpAllinds = STUDY.cluster(1,clsIdx).allinds{groupIdx,1}(icIdx);
                    tmpSetinds = STUDY.cluster(1,clsIdx).setinds{groupIdx,1}(icIdx);
                    tmpDataset = STUDY.design(STUDY.currentdesign).cell(tmpSetinds).dataset;
                    STUDY.cluster(clsIdx).topoall_group{groupIdx,1}{1,icIdx} = std_readtopo(ALLEEG, tmpDataset, tmpAllinds);
                end
            else
                for icIdx = 1:length(STUDY.cluster(1,clsIdx).setinds{1,groupIdx})
                    tmpAllinds = STUDY.cluster(1,clsIdx).allinds{1,groupIdx}(icIdx);
                    tmpSetinds = STUDY.cluster(1,clsIdx).setinds{1,groupIdx}(icIdx);
                    tmpDataset = STUDY.design(STUDY.currentdesign).cell(tmpSetinds).dataset;
                    STUDY.cluster(clsIdx).topoall_group{1,groupIdx}{1,icIdx} = std_readtopo(ALLEEG, tmpDataset, tmpAllinds);
                end
            end
        end
    end
    
% This is the case for EEGLAB15 or above.
else
    % Identify groupIdx vector.
    groupVector = zeros(length(alleegGroupLabels),1);
    for alleegIdx = 1:length(alleegGroupLabels)
        for numUniqueLabelIdx = 1:length(uniqueGroupLabels)
            currentSetIdx = find(strcmp(alleegGroupLabels,uniqueGroupLabels{numUniqueLabelIdx}));
            groupVector(currentSetIdx) = numUniqueLabelIdx;
        end
    end
    
    % Separate topos.
    for clsIdx = 2:length(STUDY.cluster)
        currentSets = STUDY.cluster(clsIdx).sets;
        for groupIdx = 1:length(uniqueGroupLabels)
            currentClusterTopoIdx = find(ismember(currentSets, find(groupVector==groupIdx)));
            currentGroupTopos = STUDY.cluster(clsIdx).topoall(1,currentClusterTopoIdx);
            if isGroupVer1orVer2 == 1
                STUDY.cluster(clsIdx).topoall_group{groupIdx,1} = currentGroupTopos;
            else
                STUDY.cluster(clsIdx).topoall_group{1,groupIdx} = currentGroupTopos;
            end
        end
    end
end
disp('Separating IC topos into groups done.')



function smoothdata = myeegfilt(data, srate, locutoff, hicutoff)
disp('Applying low pass filter (Hamming).')
disp('Transition band width is the half of the passband edge frequency.')
tmpFiltData.data   = data;
tmpFiltData.srate  = srate;
tmpFiltData.trials = 1;
tmpFiltData.event  = [];
tmpFiltData.pnts   = length(data);
tmpFiltData.nbchan = size(data,1);
filtorder = pop_firwsord('hamming', srate, hicutoff/2);
tmpFiltData_done   = pop_eegfiltnew(tmpFiltData, locutoff, hicutoff, filtorder);
smoothdata         = tmpFiltData_done.data;



function [pValue, ciList] = statPvaf(STUDY, pvafSortingIdxSelected, arglist, originalConvTopoErpAllClusters, normalizationFactors)

% (From orignal code)
frames = size(originalConvTopoErpAllClusters(:,arglist.totalPlotRange(1):arglist.totalPlotRange(end),:,:,:),2);

% (From original code) Find limits of the component selection window.
xmin = arglist.timerange(1);
xmax = arglist.timerange(2);
if any(arglist.limitcontribtime ~= 0)
    if arglist.limitcontribtime(1)<xmin
        arglist.limitcontribtime(1) = xmin;
    end
    if arglist.limitcontribtime(2)>xmax
        arglist.limitcontribtime(2) = xmax;
    end
    limframe1 = arglist.relativeFramesWithinPlotRange(1);
    limframe2 = arglist.relativeFramesWithinPlotRange(end);
    arglist.drawvertline(end+1) =  arglist.limitcontribtime(1);
    arglist.drawvertline(end+1) =  arglist.limitcontribtime(2);
else
    limframe1 = 1;
    limframe2 = frames;
end;
dataInterval = [limframe1, limframe2];

% Prepare parameters for bootstrap test.
alpha = arglist.statPvaf(1);
numIteration = arglist.statPvaf; 
designLabel = {STUDY.design(STUDY.currentdesign).variable.label};
groupFieldIdx = find(strcmp(designLabel, 'group')|strcmp(designLabel, 'session'));



%%%%%%%%%%%%%
%%% A - B %%%
%%%%%%%%%%%%%
if length(arglist.designVector)==4
    
    inputClusterIdx = arglist.clust_grandERP(pvafSortingIdxSelected);
    clustLabel = ['Cluster ' num2str(inputClusterIdx)];
    
    % Check if the difference is within or between subjects.
    if isempty(groupFieldIdx)
        numIC1 = length(STUDY.cluster(1,inputClusterIdx).topoall);
        numIC2 = length(STUDY.cluster(1,inputClusterIdx).topoall);
        topo10 = reshape(cell2mat(STUDY.cluster(1,inputClusterIdx).topoall), [67 67 numIC1]);
        topo20 = reshape(cell2mat(STUDY.cluster(1,inputClusterIdx).topoall), [67 67 numIC2]);
    else
        numIC1 = length(STUDY.cluster(1,inputClusterIdx).topoall_group{1, arglist.designVector(1,2)});
        numIC2 = length(STUDY.cluster(1,inputClusterIdx).topoall_group{1, arglist.designVector(1,4)});
        topo10 = reshape(cell2mat(STUDY.cluster(1,inputClusterIdx).topoall_group{1,arglist.designVector(1,2)}), [67 67 numIC1]);
        topo20 = reshape(cell2mat(STUDY.cluster(1,inputClusterIdx).topoall_group{1,arglist.designVector(1,4)}), [67 67 numIC2]);
    end
    
    % Combine ERPs across conditions
    ERP1 = STUDY.cluster(1,inputClusterIdx).erpdata{arglist.designVector(1,1),arglist.designVector(1,2)}(arglist.totalPlotRange,:);
    ERP2 = STUDY.cluster(1,inputClusterIdx).erpdata{arglist.designVector(1,3),arglist.designVector(1,4)}(arglist.totalPlotRange,:);
    ERP1 = ERP1(dataInterval(1):dataInterval(2),:);
    ERP2 = ERP2(dataInterval(1):dataInterval(2),:);
    combinedERP = [ERP1 ERP2];
    
    % Combine vectorized scalp topography
    topo11 = topo10(1:3:67,1:3:67,:);
    topo12 = reshape(topo11, [23^2 numIC1]);
    topo13 = topo12(~isnan(topo12));
    topo14 = reshape(topo13, [length(topo13)/numIC1 numIC1]);
    topo21 = topo20(1:3:67,1:3:67,:);
    topo22 = reshape(topo21, [23^2 numIC2]);
    topo23 = topo22(~isnan(topo22));
    topo24 = reshape(topo23, [length(topo23)/numIC2 numIC2]);
    combinedTopos = [topo14 topo24];
    
    % Set up bootstrap iteration
    selectionNumber = floor((size(topo14,2)+size(topo24,2))/2);
    if mod(selectionNumber,2)==1
        selectionNumber = selectionNumber + 1;
    end
    
    % Calculate truePvaf
    convStack1 = originalConvTopoErpAllClusters(:,arglist.pvafFrames,:,arglist.designVector(1,1),arglist.designVector(1,2));
    convStack2 = originalConvTopoErpAllClusters(:,arglist.pvafFrames,:,arglist.designVector(1,3),arglist.designVector(1,4));
    
    convSelectedCluster1 = convStack1(:,:,pvafSortingIdxSelected);
    convSelectedCluster2 = convStack2(:,:,pvafSortingIdxSelected);
    
    allClusterEnvtopo           = sum((convStack1 - convStack2),3);
    allClusterEnvtopoVar        = mean(var(allClusterEnvtopo));
    clusterSubtractedEnvtopo    = allClusterEnvtopo-((convSelectedCluster1-convSelectedCluster2));
    clusterSubtractedEnvtopoVar = mean(var(clusterSubtractedEnvtopo));
    trueDeltaPvaf = 100-100*(clusterSubtractedEnvtopoVar/allClusterEnvtopoVar);
    
    % Compute surrogate pvaf
    surroStack1 = convStack1;
    surroStack2 = convStack2;
    
    % Enter the main loop
    surroDeltaPVAF = zeros(numIteration,1);
    timeKeeper     = tic;
    timeLapseIdx   = 0;
    for iterationIdx = 1:numIteration
        
        if mod(iterationIdx,100)==0
            timeLapseIdx = timeLapseIdx+1;
            timeLapse(timeLapseIdx) = toc(timeKeeper);
            disp(sprintf('%.0f/%.0f, %.1f s remaining...', iterationIdx, numIteration, median(timeLapse)*(numIteration-iterationIdx)/100))
            timeKeeper = tic;
        end
        
        % Divide the bootstrap index equiprobably between conditions (by Zhou, Gao, Hui, 1997)
        bootIdxIdx = randi(size(combinedERP,2), size(combinedERP,2), 1);
        bootIdx1 = bootIdxIdx(1:size(ERP1,2));
        bootIdx2 = bootIdxIdx(size(ERP1,2)+1:end);
        
        % Compute surrogate ERP
        surroERP1 = combinedERP(:,bootIdx1);
        surroERP2 = combinedERP(:,bootIdx2);
        
        % Compute surrogate Topos
        surroTopo1 = combinedTopos(:,bootIdx1);
        surroTopo2 = combinedTopos(:,bootIdx2);
        
        % Convolve topo x ERP, then normalize projection
        surroConv1 = surroTopo1*surroERP1'...
            /normalizationFactors(1,1,pvafSortingIdxSelected,1, arglist.designVector(1,2));
        surroConv2 = surroTopo2*surroERP2'...
            /normalizationFactors(1,1,pvafSortingIdxSelected,1, arglist.designVector(1,4));
        
        % Stack surrogate topo_x_ERP
        surroStack1(:,:,pvafSortingIdxSelected) = surroConv1;
        surroStack2(:,:,pvafSortingIdxSelected) = surroConv2;
        surroStackSummed1 = sum(surroStack1,3);
        surroStackSummed2 = sum(surroStack2,3);
        
        % Compute surrogate PVAF.
        surroAllClusterEnvtopo           = sum((surroStackSummed1 - surroStackSummed2),3);
        surroAllClusterEnvtopoVar        = mean(var(surroAllClusterEnvtopo));
        surroClusterSubtractedEnvtopo    = allClusterEnvtopo-((surroConv1-surroConv2));
        surroClusterSubtractedEnvtopoVar = mean(var(surroClusterSubtractedEnvtopo));
        surroDeltaPVAF(iterationIdx)     = 100-100*(surroClusterSubtractedEnvtopoVar/surroAllClusterEnvtopoVar);
    end
    pValue  = stat_surrogate_pvals(surroDeltaPVAF, trueDeltaPvaf, 'upper');
    confInt = stat_surrogate_ci(surroDeltaPVAF, 0.05, 'upper')';
        

    
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (A - B) - (C - D) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
elseif length(arglist.designVector)==8
    
    inputClusterIdx = arglist.clust_grandERP(pvafSortingIdxSelected);
    clustLabel = ['Cluster ' num2str(inputClusterIdx)];
    
    % Check if the difference is within or between subjects.
    if isempty(groupFieldIdx)
        numIC1 = length(STUDY.cluster(1,inputClusterIdx).topoall);
        numIC2 = length(STUDY.cluster(1,inputClusterIdx).topoall);
        topo10 = reshape(cell2mat(STUDY.cluster(1,inputClusterIdx).topoall), [67 67 numIC1]);
        topo20 = reshape(cell2mat(STUDY.cluster(1,inputClusterIdx).topoall), [67 67 numIC2]);
    else
        numIC1 = length(STUDY.cluster(1,inputClusterIdx).topoall_group{1, arglist.designVector(1,2)});
        numIC2 = length(STUDY.cluster(1,inputClusterIdx).topoall_group{1, arglist.designVector(1,4)});
        topo10 = reshape(cell2mat(STUDY.cluster(1,inputClusterIdx).topoall_group{1,arglist.designVector(1,2)}), [67 67 numIC1]);
        if numIC2 == numIC1 % ERP1 and ERP2 from same group
            numIC2 = length(STUDY.cluster(1,inputClusterIdx).topoall_group{1, arglist.designVector(1,6)});
            topo20 = reshape(cell2mat(STUDY.cluster(1,inputClusterIdx).topoall_group{1,arglist.designVector(1,6)}), [67 67 numIC2]);
            similarGroupFlag = 1;
        else
            topo20 = reshape(cell2mat(STUDY.cluster(1,inputClusterIdx).topoall_group{1,arglist.designVector(1,4)}), [67 67 numIC2]);
            similarGroupFlag = 0;
        end
    end
    
    % combine ERP
    ERP1 = STUDY.cluster(1,inputClusterIdx).erpdata{arglist.designVector(1,1),arglist.designVector(1,2)}(arglist.pvafFrames,:);
    ERP2 = STUDY.cluster(1,inputClusterIdx).erpdata{arglist.designVector(1,3),arglist.designVector(1,4)}(arglist.pvafFrames,:);
    ERP3 = STUDY.cluster(1,inputClusterIdx).erpdata{arglist.designVector(1,5),arglist.designVector(1,6)}(arglist.pvafFrames,:);
    ERP4 = STUDY.cluster(1,inputClusterIdx).erpdata{arglist.designVector(1,7),arglist.designVector(1,8)}(arglist.pvafFrames,:);
    
    combinedERP1 = [ERP1 ERP2];
    combinedERP2 = [ERP3 ERP4];
    
    clear ERP*
    
    % combine topo
    topo11 = topo10(1:3:67,1:3:67,:);
    topo12 = reshape(topo11, [23^2 numIC1]);
    topo13 = topo12(~isnan(topo12));
    topo14 = reshape(topo13, [length(topo13)/numIC1 numIC1]);
    
    topo21 = topo20(1:3:67,1:3:67,:);
    topo22 = reshape(topo21, [23^2 numIC2]);
    topo23 = topo22(~isnan(topo22));
    topo24 = reshape(topo23, [length(topo23)/numIC2 numIC2]);
    
    if similarGroupFlag == 0
        combinedTopo1 = [topo14 topo24];
        combinedTopo2 = [topo14 topo24];
    elseif similarGroupFlag == 1
        combinedTopo1 = [topo14 topo14];
        combinedTopo2 = [topo24 topo24];
    end
    
    clear topo1* topo2*
    
    % calculate trueDeltaPvaf
    convStack1 = originalConvTopoErpAllClusters(:,arglist.pvafFrames,:,arglist.designVector(1,1),arglist.designVector(1,2));
    convStack2 = originalConvTopoErpAllClusters(:,arglist.pvafFrames,:,arglist.designVector(1,3),arglist.designVector(1,4));
    convStack3 = originalConvTopoErpAllClusters(:,arglist.pvafFrames,:,arglist.designVector(1,5),arglist.designVector(1,6));
    convStack4 = originalConvTopoErpAllClusters(:,arglist.pvafFrames,:,arglist.designVector(1,7),arglist.designVector(1,8));
    
    convSelectedCluster1 = convStack1(:,:,pvafSortingIdxSelected);
    convSelectedCluster2 = convStack2(:,:,pvafSortingIdxSelected);
    convSelectedCluster3 = convStack3(:,:,pvafSortingIdxSelected);
    convSelectedCluster4 = convStack4(:,:,pvafSortingIdxSelected);

    allClusterEnvtopo           = sum((convStack1 - convStack2) - (convStack3 - convStack4),3);
    allClusterEnvtopoVar        = mean(var(allClusterEnvtopo));
    clusterSubtractedEnvtopo    = allClusterEnvtopo-((convSelectedCluster1-convSelectedCluster2)-(convSelectedCluster3-convSelectedCluster4));
    clusterSubtractedEnvtopoVar = mean(var(clusterSubtractedEnvtopo));
    trueDeltaPvaf = 100-100*(clusterSubtractedEnvtopoVar/allClusterEnvtopoVar);
    
    % Prepare bootstrap. Initialize envtopo stack data and loop
    surroStack1 = convStack1;
    surroStack2 = convStack2;
    surroStack3 = convStack3;
    surroStack4 = convStack4;
    
    clear conv*
    
    surroDeltaPVAF = zeros(numIteration,1);
    timeKeeper   = tic;
    timeLapseIdx = 0;
    for iterationIdx = 1:numIteration
        
        if mod(iterationIdx,100)==0
            timeLapseIdx = timeLapseIdx +1;
            timeLapse(timeLapseIdx) = toc(timeKeeper);
            disp(sprintf('%.0f/%.0f, %.1f s remaining...', iterationIdx, numIteration, median(timeLapse)*(numIteration-iterationIdx)/100))
            timeKeeper = tic;
        end
        
        bootstrpIdx1 = randi(size(combinedERP1,2), [1 size(combinedERP1,2)]);
        bootstrpIdx2 = randi(size(combinedERP2,2), [1 size(combinedERP2,2)]);
        
        % surrogate ERP
        surroERP1 = combinedERP1(:,bootstrpIdx1(1:numIC1));
        surroERP2 = combinedERP1(:,bootstrpIdx1(numIC1+1:end));
        surroERP3 = combinedERP2(:,bootstrpIdx2(1:numIC2));
        surroERP4 = combinedERP2(:,bootstrpIdx2(numIC2+1:end));
        
        % surrogate Topos
        surroTopo1 = combinedTopo1(:,bootstrpIdx1(1:numIC1));
        surroTopo2 = combinedTopo1(:,bootstrpIdx1(numIC1+1:end));
        surroTopo3 = combinedTopo2(:,bootstrpIdx2(1:numIC2));
        surroTopo4 = combinedTopo2(:,bootstrpIdx2(numIC2+1:end));
        
        % convolve topo x ERP, normalize projection
        surroConv1 = surroTopo1*surroERP1'...
            /normalizationFactors(1,1,pvafSortingIdxSelected,1, arglist.designVector(1,2));
        surroConv2 = surroTopo2*surroERP2'...
            /normalizationFactors(1,1,pvafSortingIdxSelected,1, arglist.designVector(1,4));
        surroConv3 = surroTopo3*surroERP3'...
            /normalizationFactors(1,1,pvafSortingIdxSelected,1, arglist.designVector(1,6));
        surroConv4 = surroTopo4*surroERP4'...
            /normalizationFactors(1,1,pvafSortingIdxSelected,1, arglist.designVector(1,8));
        
        % replace surrogate Topo
        surroStack1(:,:,pvafSortingIdxSelected) = surroConv1;
        surroStack2(:,:,pvafSortingIdxSelected) = surroConv2;
        surroStack3(:,:,pvafSortingIdxSelected) = surroConv3;
        surroStack4(:,:,pvafSortingIdxSelected) = surroConv4;
        
        surroStackSummed1 = sum(surroStack1,3);
        surroStackSummed2 = sum(surroStack2,3);
        surroStackSummed3 = sum(surroStack3,3);
        surroStackSummed4 = sum(surroStack4,3);
        
        % Compute surrogate PVAF.
        surroAllClusterEnvtopo           = sum((surroStackSummed1 - surroStackSummed2) - (surroStackSummed3 - surroStackSummed4),3);
        surroAllClusterEnvtopoVar        = mean(var(surroAllClusterEnvtopo));
        surroClusterSubtractedEnvtopo    = allClusterEnvtopo-((surroConv1-surroConv2)-(surroConv3-surroConv4));
        surroClusterSubtractedEnvtopoVar = mean(var(surroClusterSubtractedEnvtopo));
        surroDeltaPVAF(iterationIdx)     = 100-100*(surroClusterSubtractedEnvtopoVar/surroAllClusterEnvtopoVar);
    end
    pValue  = stat_surrogate_pvals(surroDeltaPVAF, trueDeltaPvaf, 'upper');
    confInt = stat_surrogate_ci(surroDeltaPVAF, 0.05, 'upper')';
end

% Display results
disp(sprintf('\n\n%s: p = %.4f, Delta pvaf %.1f (95%% C.I. %.1f)\n\n',...
    clustLabel, pValue, trueDeltaPvaf, confInt(2)))



% Helper function to concatenate multiple cell inputs. 08/13/2018 Makoto.
function outString = summarizeMultipleLabels(inputCells)

% Do nothing if empty input.
if isempty(inputCells)
    outString = '';
    return
    
% Do nothing if non-cell input.
elseif ~iscell(inputCells)
    outString = num2str(inputCells);
    return
end

% Concatenate multiple strings in cells.
concatenatedString = '[';
for cellIdx = 1:length(inputCells)
    concatenatedString = [concatenatedString num2str(inputCells{cellIdx})]; % In case the condition names are numeric.
    if cellIdx ~= length(inputCells)
        concatenatedString = [concatenatedString  ' & '];
    else
        concatenatedString = [concatenatedString  ']'];
    end
end
outString = concatenatedString;