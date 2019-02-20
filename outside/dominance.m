function [GenDom,AveDom,Rsqr,CoordSave] = dominance(Y,X,filename)
%   The dominance function was writen August 2010 by Stephen Broomell,
%   Florian Lorenz, and Nathaniel E. Helwig

%   Purpose: dominane_fast.m computes the partition of R^2 across a set of
%            predictors in accordance with the Azen and Budescu's (2003) 
%            Psychological Methods paper

%   Inputs:
%       1) "Y" is an [N × 1] vector containing the N observations' scores
%           on the dependent variable
%       2) "X" is an [N × K] matrix where column k (k=1:K) contains the N
%           observations' scores on the k-th independent variable
%       3) "filename" is an OPTIONAL input that specifies a filename (as a 
%           string input) which is used to create a Microsoft Excel file
%           from the output variables
%          Note: Writing Excel files can take a few extra minutes. You can
%          omit the third input to skip the Excel file writing.

%   Outputs:
%       1) "GenDom" is a [1 × K] vector providing the general dominance
%           (i.e., the R^2 partition) for the K columns of "X"
%       2) "AveDom" is a [(K+1) × (K+1)] matrix such that the first column
%           contains the row labels (0:K) that correspond to the # of
%           possible variables included in the model (note: first row and
%           0 label correspond to null model including only the intercept),
%           and the remaining K columns contain the average dominance for 
%           the K variables (averaged separately within each # of possible
%           variables included in the model)
%       3) "Rsqr" is a [M × (K+1)] matrix where the first column contains
%           the R^2 for each of the M possible subset models, and the
%           remaining K columns contain the change in R^2 associated with 
%           the addition of the k-th variable (k=1:K) to the subset model
%       4) "CoordSave" is a [M × K] matrix where each row defines the
%           variables contained in the m-th (m=1:M) subset model
%          Note: the command [CoordSave Rsqr] will combine the 3rd and 4th
%          outputs into one matrix (making it possible to readily identify
%          the variables involved in each row of "Rsqr")

%   There is no bootstrap component to this program. 

%   dominance_fast.m is currently set to allow a max of 13 predictors. 

% Track Time
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check Inputs
if nargin < 2
    error('myApp:argChk', 'Wrong number of input arguments')
elseif nargin > 3
    error('myApp:argChk', 'Wrong number of input arguments')
end
[RY,CY] = size(Y);
[RX,CX] = size(X);                    % CX is the number of variables in X

if RY ~= RX                           % Check that vectors are same length
    error('myApp:dimChk', 'Dimensions of Y and X do not match')
end
if CY ~= 1
    error('myApp:dimChk', 'Y can have only 1 column')   
end
if CX > 13                            % Max predictors set at 13 currently
    error('myApp:dimChk', 'X can have no more than 13 columns')   
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pre-allocate matrix "Rsqr"
Rindex=diag(fliplr(pascal(CX+1)));    % Vector of 'CX choose 0:CX'
CTotalR=cumsum(Rindex);               % Cumulative sum of above vector
TotalR=CTotalR(length(CTotalR));      % Total sum of above vector
Rsqr=zeros(TotalR,CX+1);              % Initialize "Rsqr" matrix


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute R^2 for each individual predictor
Rsqr(1,2:CX+1)=corr(X,Y).^2;          % Vector of R^2 for each predictor
Rsqr(2:CX+1,1)=Rsqr(1,2:CX+1)';       % Fills-in null models w/ 1 predictor


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute R^2 for all combinations of predictors (i.e., subset models)
CoordSave=zeros(TotalR,CX);           % Initialize "CoordSave" matrix
CoordSave(2:CX+1,1)=1:CX;             % Fills-in null models w/ 1 predictor
yy=Y'*Y;                              % Inner product of "Y" (used below)
ssto=yy-Y'*ones(RX)*Y./RX;            % Total sum-of-squares for the models
for j=2:CX
    Coord=combilight(CX,j);           % Calls combilight.m function (below)
    [R1,C1]=size(Coord);              % # of subset models for CX choose j
    xy=[(ones(R1,RX)*Y)'; reshape(X(:,Coord)'*Y,R1,C1)'];
        % Above line creates "xy", which is a [(C1+1) × R1] matrix where
        % column r (r=1:R1) contains Xr'*Y, where Xr represents the design 
        % matrix for the subset model defined by row r (r=1:R1) of "Coord"
    for r=1:R1
        x=[ones(RX,1) X(:,Coord(r,:))];
            % Above line creates "x", which is the [RX × (C1+1)] design 
            % matrix for the subset model defined by row r of "Coord"
        Rsqr(CTotalR(j)+r,1)=1-(yy-(((x'*x)^-1)*xy(:,r))'*xy(:,r))/ssto;
            % Above line fills-in the appropriate rows of the first column 
            % of "Rsqr" with the R^2 for the corresponding subset model
    end
    CoordSave(CTotalR(j)+1:CTotalR(j+1),:) = [Coord zeros(R1,CX-C1)];
            % Above line saves the "Coord" matrix after filling-in the
            % remaining cells with zeros
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute change in R^2 associated with addition of variables
for k=2:TotalR
    cvec=1:CX;                        % List of all variables
    cs=CoordSave(k,:);                % Coords defining subset model k
    nz=cs(cs~=0);                     % Saves the non-zero coords
    nzl=length(nz)+1;                 % # of non-zero coords plus one
    cvec(nz)=[];                      % Remove coords of vars in null model
    for m=cvec
        cvek=[sort([nz m]) zeros(1,CX-nzl)];
            % Above line creates a new coord row vector defining the full
            % model that the null model will be compared to
        chck=CoordSave==ones(TotalR,1)*cvek;
            % Above line checks to find the row of "CoordSave" that
            % corresponds to the full model that we are interested in
        Rsqr(k,m+1)=Rsqr(sum(chck,2)==CX,1)-Rsqr(k,1);
            % Above line fills-in "Rsqr" with the change in R^2 between
            % the full and null models
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create average and general dominance tables
AveDom = zeros(CX+1,CX+1);            % Initialize "AveDom" table
AveDom(1,:)=Rsqr(1,:);                % Fills-in first row of "AveDom"
AveDom(2:CX+1,1)=(1:CX)';             % Fills-in the tables' row labels
for n = 2:CX
    rs=Rsqr(CTotalR(n-1)+1:CTotalR(n),2:CX+1);
            % Above line grabs the needed chunck of "Rsqr"
    AveDom(n,2:CX+1)=sum(rs)./sum(rs~=0);
            % Above line calculate the average dominance for each variable
end
GenDom=mean(AveDom(1:CX,2:CX+1));     % Calculates the general dominance
AveDom(CX+1,2:CX+1)=GenDom;           % Adds general dominance to table

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Write to Excel File
if nargin == 3
    %Create labels X1 - xk
    Label = cell(1,CX);
    for i = 1:CX
        xtemp = ['x' num2str(i)];
        Label{1,i} = char(xtemp);
    end
    %Turn of warning that matlab is writing a new sheet
    warning off MATLAB:xlswrite:AddSheet
    
    %Write General Dominance table with labels
    xlswrite(filename, Label, 'General Dominance', 'B1')
    xlswrite(filename, {'General Dominance'}, 'General Dominance', 'A2')
    xlswrite(filename, GenDom, 'General Dominance', 'B2')
    
    %Write Average Dominance table with labels
    xlswrite(filename, Label, 'Average Dominance', 'B1')
    xlswrite(filename, AveDom, 'Average Dominance', 'A2')
    
    %Write All Dominance R^2 table with labels
    xlswrite(filename, [{'Null Model Set'}, cell(1,CX-1) , {'R^2'}, Label], 'All Dominance', 'A1')
    xlswrite(filename, [CoordSave Rsqr], 'All Dominance', 'A2')  
end

toc
end %of main program


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Subfunction Called in Line 53 (or thereabouts) Follows
function [D] = combilight(N,K)
%COMBINATOR.M  was written by Matt Fig (popkenai@yahoo.com) and can compute
%combinations and permutations, both with and without repetitions.
%COMBILIGHT is a drastically stripped-down version useful for only the
%specific case needed in this program (combinations without repetitions).
%
%Florian Lorenz, July 2010

M = double(N);  % Single will give us trouble on indexing.

if K == 1
    D =(1:N).';  % These are simple cases.
    return
elseif K == N
    D = (1:N);
    return
elseif K == 2 && N > 2  % This is an easy case to do quickly.
    BC = (M-1)*M / 2;
    F = ((M-1):-1:2);
    %subfunction CUMSUM2 from original COMBINATOR integrated here
    if isfloat(F)
        F = cumsum(F);  % For single and double, use built-in.
    else
        try
            F = cumsumall(F);  % User has the MEX-File ready?
        catch
            warning('Cumsumming by loop.  MEX cumsumall.cpp for speed.') 
            for ii = 2:size(F,1)
                F(ii,:) = F(ii,:) + F(ii-1,:); % User likes it slow.
            end
        end
    end
    id=F+1;
    D = zeros(BC,2,class(N));
    D(:,2) = 1;
    D(1,:) = [1 2];
    D(id,1) = 1;
    D(id,2) = -((N-3):-1:0);
    %subfunction CUMSUM2 from original COMBINATOR integrated here again 
    if isfloat(D)
        D = cumsum(D);  % For single and double, use built-in.
    else
        try
            D = cumsumall(D);  % User has the MEX-File ready?
        catch
            warning('Cumsumming by loop.  MEX cumsumall.cpp for speed.') 
            for ii = 2:size(D,1)
                D(ii,:) = D(ii,:) + D(ii-1,:); % User likes it slow.
            end
        end
    end
    return
end

WV = 1:K;  % Working vector.
lim = K;   % Sets the limit for working index.
inc = 1;   % Controls which element of WV is being worked on.
BC = prod(M-K+1:M) / (prod(1:K));  % Pre-allocation.
D = zeros(round(BC),K,class(N));
D(1,:) = WV;  % The first row.
for ii = 2:(BC - 1)
    if logical((inc+lim)-N) % The logical is nec. for class single(?)
        stp = inc;  % This is where the for loop below stops.
        flg = 0;  % Used for resetting inc.
    else
        stp = 1;
        flg = 1;
    end

    for jj = 1:stp
        WV(K  + jj - inc) = lim + jj;  % Faster than a vector assignment.
    end
    D(ii,:) = WV;  % Make assignment.
    inc = inc*flg + 1;  % Increment the counter.
    lim = WV(K - inc + 1 );  % lim for next run.
end
D(ii+1,:) = (N-K+1):N;           
end %of function 
