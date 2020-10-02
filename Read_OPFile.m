function OP=Read_OPFile(fn)
%Read_OPFile    Read an OP-file (TLS calibration)
%
%  OP=Read_TLS_OPFile(filename) reads the object point file whose file name
%  including path is given as filename, and stores the content of the file
%  in the structure OP with the following fields:
%
%    ID          ... numeric point id (integer >0, but treated as double
%                    for simplicity) 
%    XYZ         ... 3xn array of OP coordiates, 1 column per point,
%                    rows in order X, Y, Z, values in m
%    VCM         ... sparse 3nx3n variance-covariance matrix of XYZ(:), 
%                    i.e., of the coordinates in the sequence X1, Y1, Z1,
%                    X2, ... , Zn, values in m^2
%    sysname     ... name of coordinate system
%    hand        ... 'right' or 'left', indicating handedness of the XYZ
%                    coordinate system
%    description ... brief description as given in the file (e.g.
%                    the source of the coordinates)
%    filename    ... name of the file from which the OPs were read
%                    
%  If the filename is [] or not passed, a function will be called, which
%  opens a dialogue for file selection; OP = [] will be returned, if
%  the dialogue is cancelled. The function terminates with an error
%  message if the filename does not exist, if the file cannot be opened
%  or if it cannot be read properly.
%
% See also DemoOP.txt

% AW, 2015-05-19

% if input argument is missing, initialize as default
if nargin < 1
    fn = [];
end

% initialize output
OP = [];

% check input argument, and if necessary call file selection dialogue
if isempty(fn) || ischar(fn)
    fn = getfname(fn,1,{'*.txt','OP files (*.txt)';...
                               '*.*','All Files (*.*)' },...
                               'Select OP file',[mfilename,'_in'],1);
                           
    if isempty(fn)  % terminate the function if the user has cancelled
        return;
    end
else
    error 'Input argument must be empty or valid file name.';
end

% open the input file in text mode for reading only, read row by row,
% ignore leading/trailing whitespace, and comments starting with %, 
% parse the rest
hf = fopen(fn,'rt');    % open the file and retrieve its handle
iLineCount = 0;         % initialize the line counter (number of lines read)

OP = struct(...         % initialize the output structure
     'ID',[], ...
     'XYZ', [], ...
     'VCM', [], ...
     'sysname', [], ...
     'hand',[], ...
     'description', [], ...
     'filename', fn ... 
    );

bReadingVals = 0;       % flag: currently reading values (1 after OPSTRT
                        %       and until OP_END 
bReadingVCM  = 0;       % flag: same for vcm
bVCMstrt     = 0;       % flag: has VCSTRT been encountered before?

blksiz = 100;           % increase temporarily needed point buffers
                        % by steps of blksiz points (for performance)
npts   = 0;             % number of points read/stored so far
nbuf   = 0;             % length of buffers (nbuf-npts is the number
                        % of still available rows in the buffer; the
                        % buffer will be increased as needed)
                        
bptid = [];             % buffer which will hold the encountered point ids    
bxyz  = [];             % same for coordinates
bvcm1 = [];             % same for vectorized vcm of the points

nbufcov = 0;            % number of columns in buffer of covariances
ncovprs = 0;            % number of covariance blocks actually used

bcovid = [];            % covariance point index buffer
bcovvl = [];            % covariance buffer

while ~feof(hf)         % a loop over all the file contents
    
    hs = fgetl(hf);     % get one line of text from the file
    iLineCount = iLineCount + 1;    % increment the line counter
    
    hi = strfind(hs,'%');   % find any occurrence of % in the row
    if ~isempty(hi)     % and truncate the string at the % sign (i.e. ignore comments)
        hs(hi:end)=[];
    end
    
    hs = strtrim(hs);   % remove leading and trailing white space
    if isempty(hs)
        continue;
    end
        
    if bReadingVals     % are we reading values?
        if length(hs) >= 6 && isequal(upper(hs(1:6)),'OP_END')
            bReadingVals = 0;
        else
            % parse elements
            hi = strfind(hs,'%');
            if ~isempty(hi)
                hs(hi:end)=[];
            end

            hhid = [];
            hxyz = [];
            hsv  = [];
            jj = 0;
            while ~isempty(hs)
                
                [hh,hs]=strtok(hs,[' ',8]);
                hh=strtrim(hh);
                if isempty(hh)
                   continue; 
                end
                jj = jj + 1;
                
                switch jj
                    case 1  % the point id
                        hhid = str2double(hh);
                    case {2, 3, 4}    % the coordinates
                        hxyz(jj-1)= str2double(hh);
                    otherwise
                        hsv(jj-4)= str2double(hh);
                end
            end
            
            if ~isempty(hhid)
               if isnan(hhid) || hhid <= 0 || rem(hhid,1)
                   fclose(hf);
                   error('Line %i: point ID must be positive integer',iLineCount);
               end
               if length(hxyz) ~= 3 || any(isnan(hxyz)) || any(isinf(hxyz))
                   fclose(hf);
                   error('Line %i: coordinates of point %i invalid (expecting 3 real values)',iLineCount, hhid);
               end
               switch length(hsv)
                  case 0    % no value passed == NaN
                     hvcm = diag(NaN(3,1));

                  case 1    % standard deviation (equal for all three coords)
                     hvcm = diag(hsv*hsv*ones(3,1));
  
                  case 3    % standard deviations of X, Y, and Z
                     hvcm = diag(hsv.*hsv);
 
                  case 6    % lower triangle of vcm
                     hvcm = hsv([1 2  3; 2 4  5; 3 5  6]);

                  otherwise
                     fclose(hf);
                     error('Line %i: variance specification of point %i invalid (expecting 0, 1, 3 or 6 real values)',iLineCount, hhid);
                end

                % store the point in the remporary repository
                npts = npts + 1;
                if npts > nbuf  % need to increase the buffers
                   nbuf = nbuf+blksiz;
                   bptid(nbuf,1)=0;
                   bxyz(1:3,nbuf)=zeros(3,1);
                   bvcm1(1:9,nbuf)=zeros(9,1);
                end
                
                if ~isempty(find(bptid==hhid,1))
                    fclose(hf);
                    error('Line %i: point id %i is not unique (has already been use before in this file)',iLineCount, hhid);
                end

                bptid(npts,1) = hhid;
                bxyz(1:3,npts) = hxyz(:);
                bvcm1(1:9,npts) = hvcm(:);

            end            % a point was actually found in this line
       end              % ready reding/parsing a point's line   
        
        
    elseif bReadingVCM      % are we reading vcm values?
        if length(hs) >= 6 && isequal(upper(hs(1:6)),'VC_END')
            bReadingVCM = 0;
        else
            % parse elements
            hi = strfind(hs,'%');
            if ~isempty(hi)
                hs(hi:end)=[];
            end
            
            hhid1 = [];
            hhid2 = [];
            hhv   = [];
            
            jj = 0;
            while ~isempty(hs)
                [hh,hs]=strtok(hs,[' ',8]);
                hh=strtrim(hh);
                if isempty(hh)
                   continue; 
                end
                jj = jj + 1;
                
                if jj > 11
                    fclose(hf);
                    error('Line %i: too many covariances for points %i-%i (expecting 9 values).',iLineCount, hhid1, hhid2);
                end
                
                switch jj
                    case 1  % the first point id
                        hhid1 = str2double(hh);
                    case 2  % the seconde point id
                        hhid2 = str2double(hh);
                    otherwise
                        hhv(jj-2)= str2double(hh);
                end
            end
            
            if length(hhv)~= 9
                fclose(hf);
                error('Line %i: invalid number of covariance elements for points %i-%i (expecting 9 values). ',iLineCount, hhid1, hhid2);
            end
            
            if nbufcov < ncovprs+1    % if the buffer is not large enugh, 
                                    % increase it
               nbufcov = nbufcov+blksiz;
               bcovid(nbufcov,1:2)=zeros(1,2);
               bcovvl(1:9,nbufcov)=zeros(9,1);
            end
            
            % find the point ids in the list of point ids of the
            % coordinates
            if hhid1 <= 0 || hhid2 <= 0 || ...
                    isempty(find(bptid==hhid1,1)) || isempty(find(bptid==hhid2,1))
                fclose(hf);
                error('Line %i: invalid point ids (expecting indices of points defined in earlier in the file)',iLineCount);
            end
            if hhid1 == hhid2
                fclose(hf);
                error('Line %i: vcm block of ONE point (i.e. equal indices) must be defined in coordinate list earlier in the file.', iLineCount);
            end
            
            if ~isempty(find(bcovid(:,1)==hhid1 & bcovid(:,2)==hhid2,1)) ...
                    || ~isempty(find(bcovid(:,1)==hhid2 & bcovid(:,1)==hhid2,1))
                fclose(hf);
                error('Line %i: covariance for points %i-%i or %i-%i has alredy been defined before.',iLineCount, hhid1, hhid2, hhid2, hhid1);
            end
            
            ncovprs = ncovprs+1;
            bcovid(ncovprs,1:2)=[hhid1 hhid2];
            bcovvl(1:9, ncovprs)=hhv(:);
            
            
        end
        
        
    else
        % we expect a keyword (6 characters) and possibly more text
        if length(hs) < 6
            fclose(hf);
            error(['Error reading line %i: non-empty line', ...
                ' but not starting with valid token (line too short).'], iLineCount);
        end
        
        % find out which keyword we have
        switch(upper(hs(1:6)))
            
            case 'OPHAND'       % handedness
                if ~isempty(OP.hand)          % is thisthe first time? if not: error
                    fclose(hf);
                    error('Line %i: Handedness of coordinate system has already been defined before', iLineCount );
                end
                hs(1:6)=[];
                hi = strfind(hs,'=');    % find what we have after the '='
                if isempty(hi) || length(hs) == hi
                    fclose(hf);
                    error('Line %i: Handedness definition missing after OPHAND', iLineCount);
                end
                hs = lower(strtrim(hs((hi+1):end)));
                switch hs
                    case 'left'
                        OP.hand = 'left';
                    case 'right'
                        OP.hand = 'right';
                    otherwise
                        fclose(hf);
                        error('Line %i: Unknown OP handedness %s',iLineCount, hs);
                end
                
            case 'OPDSCR'
                if ~isempty(OP.description)
                    fclose(hf);
                    error('Line %i: OP description has already been given before.',iLineCount);
                end
                hs(1:6)=[];
                hi = strfind(hs,'=');    % find what we have after the '='
                if ~isempty(hi) && length(hs) > hi
                    OP.description = strtrim(hs((hi+1):end));
                else
                    OP.description = ' ';
                end
                
            case 'OPSYSN'
                if ~isempty(OP.sysname)
                    fclose(hf);
                    error('Line %i: OP system name has already been given before.',iLineCount);
                end
                hs(1:6)=[];
                hi = strfind(hs,'=');    % find what we have after the '='
                if ~isempty(hi) && length(hs) > hi
                    OP.sysname = strtrim(hs((hi+1):end));
                else
                    OP.sysname = ' ';
                end   
                
            case 'OPSTRT'
                if ~isempty(bptid)
                    fclose(hf);
                    error('Line %i: Keyword OPSTRT appearing for the 2nd time in the file.', iLineCount );
                else
                    bReadingVals = 1;
                end
                
            case 'OP_END'   % this is handled outside the switch, if it is detected herein, something is wrong
                fclose(hf);
                error('Line %i: OP_END apearing without previous OPSTRT.', iLineCount );
                
            case 'VCSTRT'
                if bVCMstrt
                    fclose(hf);
                    error('Line %i: Keyword VCSTRT appearing for the 2nd time in the file.', iLineCount );
                else
                    bReadingVCM = 1;
                    bVCMstrt = 1;
                end
                
            case 'VC_END'   % like with AP_END
                fclose(hf);
                error('Line %i: VC_END apearing without previous VCSTRT.', iLineCount );
                
            otherwise
                fclose(hf);
                error(['Error reading line %i: non-empty line', ...
                ' but not starting with valid token (%s).'], iLineCount, hs(1:6));
        end
        
    end
        
    
end   
fclose(hf);         % all file conents have now been read and loaded

% store data into output arrays
if npts < nbuf          % first remove unused elements from the buffers
   hhh = npts+1;
   bptid(hhh:end)=[];
   bxyz(:,hhh:end) = [];
   bvcm1(:,hhh:end) = [];
end
if ncovprs < nbufcov
    hhh = ncovprs + 1;
    bcovid(hhh:end,:) = [];
    bcovvl(:,hhh:end) = [];
end

OP.ID = bptid;         % the point ids can be stored as they are
OP.XYZ = bxyz;         % also the coordinates

% the vcm is setup as sparse matrix, therefore, the elements read from
% the input file are now associated with their row (hi) and column (hj)
% indices. For simplicity, this is done separately for the main diagonal
% block (hi1, hj1), for the off-diagonal blocks directly sepcified in the
% input file (hi2, hj2), and for the off-diagonal blocks corresponding
% to the given ones but lying in the other triangle of the entire vcm
% (such that the vcm is symmetric in the end).
hi1 = repmat(kron(ones(3,1),(1:3)'),npts,1) + kron( ((0:(npts-1))*3)', ones(9,1) );
hj1 = repmat(kron((1:3)', ones(3,1)),npts,1) + kron( ((0:(npts-1))*3)', ones(9,1) );
hv1 = bvcm1(:);

% build off-diagonal blocks using for loop (slow, but simpler - speed up
% later by vectorizing the assignment)
hnc = ncovprs;               % number of covariance pairs
hnp = length(OP.ID);         % number of points in point list
hid1 = zeros(hnc,1);         % map covariance point indices into rows of
hid2 = zeros(hnc,1);         %    OP.ID

for jj=1:hnc   % loop over all covariance pairs
    for kk=1:hnp     % for each of them, find the first and second point
        if bcovid(jj,1)== OP.ID(kk)
            hid1(jj) = kk;
        end
        if bcovid(jj,2)== OP.ID(kk)
            hid2(jj) = kk;
        end
        if hid1(jj)>0 && hid2(jj)>0
            break;
        end
    end
end

hi2 = repmat(kron(ones(3,1),(1:3)'),hnc,1) + kron( (hid1(:)-1)*3, ones(9,1) );
hj2 = repmat(kron((1:3)', ones(3,1)),hnc,1) + kron( (hid2(:)-1)*3, ones(9,1) );
hv2 = bcovvl(:);

hi3 = hj2;
hj3 = hi2;
hv3 = hv2;

hi = [hi1;hi2;hi3];
hj = [hj1;hj2;hj3];
hv = [hv1;hv2;hv3];

OP.VCM = sparse(hi,hj,hv);


end
