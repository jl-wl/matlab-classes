classdef Iteration
    %ITERATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Func % function handle
        NOut % number of output variables
        Vars % (cell array) variable values to be iterated
        Size % size of Vars
    end
    
    methods
        function it = Iteration(func,nout,varargin)
            it.Func = func;
            it.NOut = nout; % ask for explicit input because of anonymous functions
            it.Vars = varargin;
            it.Size = cellfun(@length,varargin);
        end
        
        function varargout = run(it,settings)
            f = it.Func;
            n = it.NOut;
            c = Iteration.varsProduct(it.Vars{:});
            if nargin<2
                settings = [];
            end
            
            qfield = @Utility.qfield;
            np = qfield(settings,'poolsize',1);
            bs = qfield(settings,'binsize',1); % size to be used by cellfun
            fn = qfield(settings,'filename','iter');
            st = qfield(settings,'savetemp',true); % save temporary data
            sc = qfield(settings,'singlecopy',[]); % keep single copy of output
            if nargout>n+1 || (nargout==n+1 && np==1)
                error('Too many requested output!')
            end
            
            ni = length(c{1}); % total number of iterations
            nb = ceil(ni/bs);
            pa = mat2cell((1:ni).',[bs*ones(1,nb-1),ni-bs*(nb-1)]); % index partition
            out = cell(nb,n);

            if np>1
                fnstr = Iteration.startLog(fn);
                parsave = @(ii,data) save([fnstr,'-',num2str(ii),'.mat'], ...
                                          'data','-v7.3');
                try
                    poolobj = parpool(np);
                catch ME
                    switch ME.identifier
                        case 'parallel:convenience:ConnectionOpen'
                            poolobj = gcp;
                        otherwise
                            rethrow(ME)
                    end
                end
                parfor (ib = 1:nb, np)
                    tStart = tic;
                    ci = cellfun(@(a)a(pa{ib}),c,'UniformOutput',false);
                    oi = cell(1,n);
                    [oi{:}] = cellfun(f,ci{:},'UniformOutput',false);
                    out(ib,:) = oi;
                    tElapsed = toc(tStart);
                    Iteration.writeLog(fnstr,ib,nb,tElapsed);
                    if st, parsave(ib,oi); end
                end
                delete(poolobj)
                Iteration.finishLog(fnstr);
                if st, delete([fnstr,'-*.mat']); end
            else
                for ib = 1:nb
                    tStart = tic;
                    ci = cellfun(@(a)a(pa{ib}),c,'UniformOutput',false);
                    [out{ib,:}] = cellfun(f,ci{:},'UniformOutput',false);
                    tElapsed = toc(tStart);
                    fprintf('%6d / %6d : %6.2f\n',ib,nb,tElapsed);
                end
            end
            
            if nargout<n+1
                varargout = arrayfun(@(a)vertcat(out{:,a}), 1:nargout, ...
                                     'UniformOutput', false);
            else % nargout==n+1 && np>1
                varargout(1:n) = arrayfun(@(a)vertcat(out{:,a}), 1:n, ...
                                          'UniformOutput', false);
                varargout{n+1} = fnstr;
            end
            
            if ~isempty(sc)
                sc = sc(sc<=min(nargout,n));
                for iout = sc
                    oi = varargout{iout};
                    varargout{iout} = oi{1};
                end
            end
        end
        
        function varargout = reshapeData(it,varargin)
        % convert data to multi-dimentional numerical array
            for ii = 1:nargin-1
                data = varargin{ii};
                s = size(data{1});
                data = reshape(cell2mat(data),[s(1),it.Size,s(2:end)]);
                if prod(s(2:end))>1
                    nit = length(it.Size);
                    data = permute(data,[1,nit+(2:length(s)),1+(1:nit)]);
                end
                varargout{ii} = data; %#ok<AGROW>
            end
        end
    end
    
    methods (Static)
        function c = varsProduct(varargin)
        % Cartesian product of variables, stored in one cell array
            l = cellfun(@length,varargin);
            N = length(l);
            v0 = cellfun(@(a)reshape(a,[],1),varargin,'UniformOutput',false);
            rep = @(a,n)repmat(repelem(a,prod(l(1:n-1)),1),prod(l(n+1:end)),1);
            c = cellfun(rep,v0,num2cell(1:N),'UniformOutput',false);
            for n = 1:N
                if isnumeric(c{n})
                    c{n} = num2cell(c{n});
                end
            end
        end
        
        function fnstr = startLog(filename)
            fnstr = [filename,Utility.ct];
            fileID = fopen([fnstr,'.log'],'w');
            fprintf(fileID,'Start at %s\n',num2str(clock,'%02.0f'));
            fclose(fileID);
        end

        function writeLog(fnstr,ii,nn,tt)
            fileID = Utility.fopenwait([fnstr,'.log'],'a');
            fprintf(fileID,'%6d / %6d : %6.2f\n',ii,nn,tt);
            fclose(fileID);
        end

        function finishLog(fnstr)
            fileID = fopen([fnstr,'.log'],'a');
            fprintf(fileID,'End at %s',num2str(clock,'%02.0f'));
            fclose(fileID);
        end
    end
end

