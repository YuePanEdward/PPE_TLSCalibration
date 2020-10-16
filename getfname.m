function fn = getfname(fn,mustexist,fspec,caption,progname,store);
%fn = getfname(fn,mustexist,fspec,caption,progname,store);
%  Determine and check a filename. If fn is empty, a dialog is invoked
%  for file selection, using the file specification fspec (see uigetfile
%  for more information) and the given dialog 'caption'.
%
%  If a program name (progname) is passed, the function tries to determine
%  which directory was last used in association with this program name, and
%  if it finds one, opens the dialog starting in that directory. The
%  association information is handled by the function 'workingdir' which
%  accesses a special ASCII file.
%
%  If the flag 'store' is passed and is non-empty and non-zero, the
%  directory from which a file is chosen is associated with the given 
%  progam name for later use.
%
%  If mustexist is non-empty and non-zero, the filename specified must
%  identify an existing file. Should the selected file not exist, the
%  function returns with a corresponding warning.
%
%  If the filename is passed as a non-empty argument (fn), it is essentially
%  passed through the function. However, the test for existence and the
%  association with a program name may be performed, according to the
%  parameters 'mustexist' and 'store'.
%
%  If the returned value is [], the user has cancelled the operation 
%  or an error occurred.

% created      : A. Wieser
% last modified: AW 3/1/2004
% non-standard functions needed: workingdir


% check/complete the input arguments
if nargin < 1
   error('Too few arguments.')
end
if nargin < 6  % initialize missing parameters as []
   store = [];
   if nargin < 5
      progname = [];
      if nargin < 4
         caption = [];
         if nargin < 3
            fspec = [];
            if nargin < 2
               mustexist = [];
            end
         end
      end
   end
end

% replace empty parameters by default values
if isempty(mustexist)
   mustexist = 0;
else
   if ~isequal(mustexist,0)
      mustexist = 1;
   end   
end
if isempty(fspec)
   fspec = {'*.*','All files (*.*)'};
end
if isempty(caption)
   caption = 'Select file';
end
if isempty(store)
   store = 0;
end


% interactively get a filename if there is none (switch to last directory,
% where a file has been taken from)
if isempty(fn)
   save_dir = pwd;                     % remember current directory
   hwd = workingdir(progname);         % get last directory associated with this program (if any...)
   if ~isempty(hwd)                    % if there is anything -> change to that directory
      cd(hwd);
   end

   if ~mustexist
      [fn,pn] = uiputfile(fspec,caption);
   else
      [fn,pn] = uigetfile(fspec,caption);
   end
                           
   if isequal(fn,0) | isequal(pn,0)    % user has cancelled
      cd(save_dir);                    % restore previous working directory
      fn = [];
      return;                          % and exit the whole function
   end
      
   fn = fullfile( pn, fn );            % otherwise compose full filename
   cd(save_dir);                       % and restore previous working directory
   
else  % i.e. fn was not empty

   % check whether the file name is a string
   if ~ischar(fn)
      error 'Filename must be a string!'
   end

   % if the file needs to exist, check for existence
   if mustexist
      % try to open it
      try
         fid = fopen(fn,'rb');
      catch
         hs = lasterror;
         error(sprintf('Cannot open file %s: %s', fn, hs.message));
      end   

      if fid < 0
          error('Cannot open file %s: %s', fn, 'File or path does not exist?');
      end    
      try
         fclose(fid);
      catch
         hs = lasterror;
         error(sprintf('Cannot access file %s: %s', fn, hs.message));
      end   
      
      
   end
   
end   


% if the path is to be stored, do this now
if ~isequal(store,0)
   if isempty(progname)
      warning('Cannot associate pathname with program: program name missing.');
   end
   pathstr = fileparts(fn);
   if isempty(pathstr)
      pathstr = pwd;
   end   
   workingdir(progname,pathstr);
end
   
   
   

