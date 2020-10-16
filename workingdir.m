function dirname = workingdir(funcname,dirname)
%dirname = workingdir(funcname)
%   Returns the working directory associated with the program
%   funcname, if there is any record of this program in
%   %MATLABHOME%\toolbox\local\awdirs.mat and if this record
%   points to a valid directory.
%
%workingdir(funcname,dirname)
%   Associates the directory dirname with the program
%   funcname (the repository is by default in %MATLABHOME%\toolbox\local\awdirs.mat).

% created      : A. Wieser
% last modified: AW 27/09/2018 (warning if path to repository does not exist)


    % what do we have to do?
    if nargin < 2                 % get previously used directory
       dirname = lfn_gdn(funcname);
       if ~isempty(dirname) && ~exist(dirname,'dir')
          dirname = [];
       end
    else                          % save directory 
       lfn_sdn(funcname,dirname);
    end
end


% **********************************************************************
% local function: get directory associated with funcname s
function y = lfn_gdn(s)

   % get function list and directory list
   [funclist,dirlist]=lfn_gfdl;
   
   % search the funclist for funcname (ignoring upper/lower case differences)
   s = lower(s);
   y = [];
   for i=1:length(funclist)
      if ischar(funclist{i}) && isequal(lower(funclist{i}),s)
         if length(dirlist) >= i
            y = dirlist{i};
            return;
         end
      end   
   end
end

% local function: store directory d associated with funcname s
function lfn_sdn(s,d)

   % get function list and directory list
   [funclist,dirlist]=lfn_gfdl;
   
   % search the funclist for funcname (ignoring upper/lower case differences)
   s = lower(s);
   ind = 0;
   nfunc = length(funclist);
   for i=1:nfunc
      if ischar(funclist{i}) && isequal(lower(funclist{i}),s)
         ind = i;
         break;
      end   
   end
   
   % insert the function name (if it is not part of funlist yet)
   if ind == 0 % i.e. the entry was not found => create it
      ind = nfunc + 1;
      funclist{ind} = lower(s);
   end
   
   % set the directory name
   dirlist{ind} = d;
   
   % and save the function list
   fn0 = lfn_gpthf;         % start with getting the path
   
   if ~exist(fn0,'dir')     % if the predefined path does not exist, warn the user
       warning('Path to store/retrieve information does not exist (''%s'').\nConsider changing in lfn_gpthf within this file.',fn0);
   end
   
   fn = [fn0,filesep,'awdirs.mat'];
   if exist(fn,'file')
      save(fn,'funclist','dirlist','-append');
   else
      save(fn,'funclist','dirlist');      
   end
end


% local function: get function and directory lists from repository
function [funclist,dirlist]=lfn_gfdl
   
   % initialize lists (in case the repository does not exist)
   funclist = [];
   dirlist = [];
      
   % try to locate the file
   fn = [lfn_gpthf,filesep,'awdirs.mat'];
   
   % if it exists
   if exist(fn,'file')
      % load all data from this file
      h = load(fn);
      
      % and extract the funclist and dirlist cell arrays therefrom
      if isfield(h,'funclist') && isfield(h,'dirlist')
         funclist = h.funclist;
         dirlist  = h.dirlist;
      end
   end
end
   
   
% local function: get path to file containing directory names (for
% performance using a global variable might be better)
function path_name=lfn_gpthf
    % by default: %MATLABHOME%\toolbox\local - if not available on your
    % computer, change here
    path_name = [matlabroot,filesep,'toolbox',filesep,'local'];
end
