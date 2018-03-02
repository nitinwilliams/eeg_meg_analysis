% InfoMap method from (Roswall and Bergstrom, 2008)
%
% Input
%   - adj: adjacency matrix
%
% Output
%   - com: communities (listed for each node)
%
% Author: Erwan Le Martelot (modified by Nitin Williams, date: 20/07/2017)
% Date: 16/06/11

function [com] = infomap(adj)

    % Set the path and command line name
    dir_path = '/usr/local/infomap/'; % should be set to 'infomap' path (also create a directory in this path called '1')
    command = 'Infomap';

    command = [dir_path,command];
    command = [command,' -u -2 --hard-partitions --clu --silent adj.net',' 1'];

    % Get edges list and output to file
    adj2pajek(adj,'adj','.');
 
    % Call community detection algorithm
    system(command);
       
    % Load data and create data structure
    A=importdata('1/adj.clu');
    DAT=A.data;
    [~,reorder]=sort(DAT(:,1),'ascend');   
    com=DAT(reorder,2);
    
end