classdef Utility
%TODO:
%   - use pagemtimes in einsum
    
    properties (Constant)
        PhysConst = struct('c', 299792458, ... % speed of light; m/s
                           'e', 1.6021766208e-19, ... % elementary charge; A.s
                           'eps0', 8.854187817620e-12, ... % vacuum permittivity; A^2.s^4/(kg.m^3)
                           'mu0', 1.256637061435917e-06, ... % vacuum permeability; kg.m/(A^2.s^2)
                           'h', 6.626070040e-34, ... % Planck constant; kg.m^2/s
                           'hbar', 1.054571800e-34, ...
                           'k', 1.38064852e-23, ... % Boltzmann constant; kg.m^2/(s^2.K)
                           'me', 9.10938356e-31, ... % electron rest mass; kg
                           'muB', 9.274009994e-24, ... % Bohr magneton; A.m^2
                           'navo', 6.022140857e23 ... % Avogadro constant; mol
                          );
        PauliMats = {[0, 1; 1, 0], ...
                     [0, -1i; 1i, 0], ...
                     [1, 0; 0, -1], ...
                     [1, 0; 0, 1]};
        DiracMats = {kron([0, 1; 1, 0],[0, 1; 1, 0]), ... % alpha1
                     kron([0, 1; 1, 0],[0, -1i; 1i, 0]), ... % alpha2
                     kron([0, 1; 1, 0],[1, 0; 0, -1]), ... % alpha3
                     kron([1, 0; 0, -1],[1, 0; 0, 1]), ... % beta
                     kron([0, -1i; 1i, 0],[1, 0; 0, 1])}; % alpha123*beta
        orbitmap = containers.Map({'s', 'p', 'd', 'f'},[0, 1, 2, 3]);
        bondmap = containers.Map({'sigma', 'pi', 'delta', 'phi'},[0, 1, 2, 3]);
    end
    
    methods (Static)
        function p = lorentzian(dx,g)
        % lorentzian (cauchy) distribution function, normalized to 1
        % g (scalar) is the FWHM
            p = g./(2*pi*(dx.^2+g^2/4));
        end
        
        function w = dfden(de,T)
        % minus derivative of Fermi-Dirac distribution 1/(1+exp(de/T)) wrt de
        % T (scalar) is the temperature
            w = 0.5./(T*(cosh(de/T)+1));
        end
        
        function C = chop(A,p)
            if nargin<2
                p = 1e-12;
            end

            if isreal(A)
                C = A.*(abs(A)>abs(p));
            else
                C = real(A).*(abs(real(A))>abs(p))+1i*imag(A).*(abs(imag(A))>abs(p));
            end
        end
        
        function result = iif(cond,t,f)
        %IIF implements a ternary operator
            if isscalar(cond) 
                if cond 
                    result = t;
                else
                    result = f;
                end
            else
              result = (cond).*t + (~cond).*f;
            end
        end
        
        function n = countperm(v,direction)
        % count permutation needed to sort vector to 'ascend'/'descend'
            if nargin<2
                direction = 'ascend'; %default
            end
            switch direction(1)
                case 'a'
                    cfun = @(ii)sum(v(1:ii-1)>v(ii));
                case 'd'
                    cfun = @(ii)sum(v(ii+1:end)>v(ii));
                otherwise
                    error('unknow direction')
            end
            n = sum(arrayfun(cfun,1:length(v)));
        end
        
        function [a, p] = parts(volumns)
        % all possible ways to partition an index set according to volumns,
        % or, put indexed balls into indexed jars of given volumns,
        % assuming total volumn equals number of balls (if not, consider a
        % virtual jar or virtual balls)
        % output: a (array), all assignments of balls to jar indices;
        %         p (cell array), all ball indices assigned to each jar
            l = length(volumns);
            a = unique(perms(repelem((1:l).',volumns)),'rows');
            f = @(ii)reshape(mod(find(a'==ii)-1,size(a,2))+1,volumns(ii),[])';
            p = arrayfun(f,1:l,'UniformOutput',false);
        end
        
        function spm = spd(v,d,pbc)
        % spd constructs a n-by-n sparse matrix with v being the d-diagonal.
        % without periodical boundary condition, v has rank (n-abs(d));
        % with pbc, v has rank n and starts with hopping FROM site 1 (if d~=0).
            if size(v,2)>1
                error('v must be a column vector!');
            end
            
            if nargin<2
                d = 0;
            end

            if nargin<3
                pbc = false;
            end

            l = length(v);
            if ~pbc
                if d>0
                    v = [zeros(abs(d),1);v];
                elseif d<0
                    v = [v;zeros(abs(d),1)];
                end
                n = l+abs(d);
                spm = spdiags(v,d,n,n);
            else
                n = l;
                spm = spdiags(v,d,n,n);
                if d~=0
                    d2 = -sign(d)*(n-abs(d));
                    spm = spm + spdiags(v,d2,n,n);
                end
            end
        end
        
        function p = mdot(A,B,func)
        % matrix-vector dot, also possible vector dot matrices
            if nargin<3
                func = @mtimes;
            end
            if isnumeric(A) && isnumeric(B)
                if nargin<3
                    p = dot(A,B);
                else
                    p = sum(arrayfun(func,A,B));
                end
                return
            elseif isnumeric(A) && iscell(B)
                p = func(A(1),B{1});
                for ii = 2:length(B)
                    p = p+func(A(ii),B{ii});
                end
            elseif iscell(A) && isnumeric(B)
                p = func(A{1},B(1));
                for ii = 2:length(A)
                    p = p+func(A{ii},B(ii));
                end
            elseif iscell(A) && iscell(B)
                p = func(A{1},B{1});
                for ii = 2:length(A)
                    p = p+func(A{ii},B{ii});
                end
            else
                error('Illegal input types!')
            end
        end
        
        function s = sumall(a,ds)
            s = a;
            if nargin<2
                while ~isscalar(s)
                    s = squeeze(sum(s));
                end
            else
                ds = sort(ds);
                for d = ds
                    s = sum(s,d);
                end
                s = squeeze(s);
            end
        end
        
        function a = kronall(varargin)
            if nargin>1
                a = 1;
                for ii = 1:nargin
                    a = kron(a,varargin{ii});
                end
            elseif nargin==1
                a = varargin{1};
            else
                error('At least one input!')
            end
        end
        
        function a = operall(op,varargin) % caution: NOT efficient
        % operations (plus, times, etc.) on a collection
            if nargin<2
                a = [];
            else
                a = varargin{1};
                for ii = 2:length(varargin)
                    a = op(a,varargin{ii});
                end
            end
        end
        
        function ct_string = ct
        % get current time string
            ct_string = strrep(num2str(clock,'%02.0f'),'  ',''); %#ok<CLOCK>
        end
        
        function fileID = fopenwait(varargin)
            fileID = -1;
            while fileID<0
                pause(0.01)
                fileID = fopen(varargin{:});
            end
        end
        
        function varargout = eigen(varargin)
        % input can be the same as eig or eigs
            m = varargin{1};
            d = size(m,1);
            large_d = 10000;

            if nargin>1
                n = varargin{2};
                if n<d
                    feig = @eigs;
                    varargin{1} = sparse(m);
                else
                    if d>large_d
%                         button = questdlg('Large matrix. Continue full?');
%                         switch button
%                             case {[],'No','Cancel'}
%                                 varargout = cell(nargout);
%                                 return
%                         end
                        warning(['Use full large matrix (dim=', num2str(d),')!'])
                    end
                    feig = @eig;
                    varargin = {full(m)};
                end
            else % when there is only one input, always use full matrix
                if d>large_d
                    warning(['Use full large matrix (dim=', num2str(d),')!'])
                end
                feig = @eig;
                varargin = {full(m)};
            end

            if nargout<2
                evals = feig(varargin{:});
                if mean(abs(imag(evals)))<1e-10
                    varargout{1} = sort(real(evals));
                else % nonhermitian
                    varargout{1} = sort(evals);
                end
            elseif nargout==2
                [V,D] = feig(varargin{:});
                evals = diag(D);
                if mean(abs(imag(evals)))<1e-10
                    [evals,indx] = sort(real(evals));
                    dim = length(evals);
                    orth_err = norm(V'*V-speye(dim),'fro')/dim;
                    if orth_err>1e-3
                        warning(['Orthogonalization Error = ',num2str(orth_err)])
                    end
                else % nonhermitian
                    [evals,indx] = sort(evals);
                end
                varargout{1} = V(:,indx);
                varargout{2} = evals;
            elseif nargout==3 && ~ishermitian(varargin{1})
                [V,D,W] = eig(varargin{:}); % A = V*D*W'
                WV = W'*V;
                WVd = diag(WV);
                WVerr = norm(WV-diag(WVd));
                if WVerr>1e-12
                    warning(['Normalization Error(W*V) = ', num2str(WVerr)])
                end
                V = V*diag(1./WVd);
                evals = diag(D);
                [evals,indx] = sort(evals);
                varargout{1} = V(:,indx);
                varargout{2} = evals;
                varargout{3} = W(:,indx);
            else
                [varargout{1:nargout}] = feig(varargin{:});
            end
        end
        
        function [l, nbs] = linInterp(l0,nb)
        % interpolate a multidimensional list (l0) with number of bins (nb)
        % inserted in each pair of adjacent slices defined wrt 1st dimension
            s = size(l0);
            n0 = s(1);

            if n0<2
                l = l0;
                return
            end

            if length(nb)<n0-1
                dl = diff(l0,1,1).';
                nb = round(2*nb*sqrt(sum(dl.^2,1))); % nb points per unit length
            %     nb = nb(1)*ones(1,n0-1);
            else
                nb = reshape(nb(1:n0-1),1,[]);
            end

            l0 = reshape(l0,n0,[]);
            l = zeros(sum(nb)+1,size(l0,2));
            offset = 0;
            for ii = 1:n0-1
                nbii = nb(ii);
                d = (l0(ii+1,:)-l0(ii,:))/nbii;
                l(offset+(1:nbii),:) = kron(ones(nbii,1),l0(ii,:))...
                                       +kron((0:nbii-1)',d);
                offset = offset+nbii;
            end
            l(offset+1,:) = l0(n0,:);
            s(1) = size(l,1);
            l = reshape(l,s);

            if nargout>1
                nbs = nb;
            end
        end
        
        function p = cartprod(varargin)
        % Cartesian product of multiple sets (arrays or cell arrays)
            ls = cellfun(@length,varargin);
            ms = arrayfun(@(ii)prod(ls(ii+1:end)),1:length(ls),'UniformOutput',false);
            ns = arrayfun(@(ii)prod(ls(1:ii-1)),1:length(ls),'UniformOutput',false);
            c = cellfun(@(a,m,n)repmat(repelem(reshape(a,[],1),m,1),n,1),...
                        varargin,ms,ns,'UniformOutput',false);
            p = cat(2,c{:});
        end
        
        function [v, l] = ndgridVecs(mins,maxs,steps)
        % generate ndgrid and put each output in one column of a matrix v
        % l is the length of vector in each grid dimension
            if nargin<3
                steps = 1;
            end
            if isscalar(steps)
                steps = repmat(steps,size(mins));
            end

            c1 = arrayfun(@(a,b,s)(a:s:b),mins,maxs,steps,'UniformOutput',false);
            l = cellfun(@length,c1);

            c2 = cell(1,length(mins));
            [c2{:}] = ndgrid(c1{:});
            v = cell2mat(cellfun(@(a)reshape(a,[],1),c2,'UniformOutput',false));
        end
        
        function val = qfield(s,field_name,default_value)
            if isfield(s,field_name)
                val = s.(field_name);
            else
                val = default_value;
            end
        end
        
        function ind = findvec(v0,vs,tol)
            if nargin<3
                tol = 1e-12;
            end
            if size(v0,1)==1 && size(v0,2)==size(vs,2)
                ind = find(sum(abs(vs-repmat(v0,size(vs,1),1)),2)<tol);
            elseif size(v0,2)==1 && size(v0,1)==size(vs,1)
                ind = find(sum(abs(vs-repmat(v0,1,size(vs,2))),1)<tol);
            else
                error('Invalid input vector to be searched!')
            end
        end
        
        function ind = findstruct(s1,s2)
            assert(isstruct(s1) && isstruct(s2), 'Inputs must be struct!')
            if isscalar(s1)
                ind = find(arrayfun(@(s)isequal(s,s1),s2));
            else
                ind = arrayfun(@(s)find(arrayfun(@(ss)isequal(ss,s),s2)),s1,...
                      'UniformOutput',false);
                if all(cellfun(@(ii)isscalar(ii),ind))
                    ind = cell2mat(ind);
                end
            end
        end
        
        function ind = findperm(A, B, tol)
        % find permutation s.t. A(ind)==B
            if nargin<3
                tol = 1e-12;
            end
            [A, iA] = sort(A);
            [B, iB] = sort(B);
            if norm(A-B)>tol
                error('A and B are not related by permutation!')
            end
            ind(iB) = iA;
        end
        
        function [sp, p] = pf(A,U)
        % To compute the sign and value of Pfaffian of a skew-symmetric matrix
        % using the Parlett-Reid algorithm.

            assert(ismatrix(A), 'argument must be a matrix')
            assert(size(A,1)==size(A,2), 'argument is not skew-symmetric')
            % make sure input is skew-symmetric
            if nargin>1
                A = U'*A*U;
            end
            assert(norm(A+A.','fro')<1e-14*size(A,1), 'argument does not seem skew-symmetric')

            N = size(A,1);

            if( mod(N,2)==1 )
                sp = 0;
                p = 0.0;
                return
            end
            
            sp = 1;
            p = 1.0;
            for k = 1:2:N-2
                % First, find the largest entry in A(k+1:N,k) and
                % permute it to A(k+1,k)
                [~, kp] = max(A(k+1:N, k),[],"ComparisonMethod","abs");
                kp = kp + k;

                % Check if we need to pivot
                if kp ~= k+1
                    % interchange rows k+1 and kp
                    A([k+1,kp],k:N) = A([kp,k+1],k:N);

                    % interchange columns k+1 and kp
                    A(k:N,[k+1,kp]) = A(k:N,[kp,k+1]);

                    % every interchange gives a "-" for the determinant of the
                    % permutation
                    sp = -sp;
                    p = -p;
                end

                sp = sp * sign(A(k,k+1));
                if nargout>1
                    p = p * A(k,k+1);
                end

                if A(k+1,k) ~= 0
                    m = (A(k+2:N,k)/A(k+1,k))*A(k+2:N,k+1).';
                    A(k+2:N,k+2:N) = A(k+2:N,k+2:N) + m - m.';
                end
            end

            sp = sp * sign(A(N-1,N));
            if nargout>1
                p = p * A(N-1,N);
            end
        end
        
        function C=einsum(A,B,iA,iB)
        % Efficiently calculates tensor contraction in any dimension. It uses
        % matlab's matrix multiplication so parallelized and optimized. The
        % Einstein summation is inspired by NumPy's syntax, but is not identical. 
        %
        % Usage: einsum(A,B,s) 
        %        einsum(A,B,iA,iB)
        %
        % Calculates the contraction of A and B which are n-dimensional tensors.
        % The contraction is specified either by s, in Einstein notation, or by two
        % vectors, iA and iB which list the indices to contract on for each tensor.
        %
        % Example:
        %     A=rand(7,4,5);
        %     B=rand(5,7);
        % To contract the 1st dimension of A with the 2nd dimension of B, use
        %   einsum(A, B, 'ijk,li->jkl') OR einsum(A, B, 1, 2)
        % The result will be of size [4,5,5]. 
        %
        % To contract the 3rd dimension of A with the 1st dimension of B, use
        %   einsum(A, B, 'ijk,kl->ijl') OR einsum(A, B, 3, 1)
        % The result will be of size [7,4,7]. 
        %
        % To do both contractions at once, use
        %   einsum(A,B,'ijk,ki->j') OR einsum(A, B, [1 3], [2 1])
        %
        % Using the iA, iB it is not possible to specify the order of dimensions
        % in the output, they're just in the same order as the inputm with the
        % contracted dimensions omitted.
        %
        % Author: Yohai Bar-Sinai 

            sA=size(A);
            sB=size(B);
            if nargin==3
                [iA, iB, final_permutation]=parse(iA, sA, sB);
            else
                final_permutation=false;
            end

            if length(iA)~=length(unique(iA)) || length(iB)~=length(unique(iB))
                error('each dimension should appear only once.')
            end

            if ~isequal(sA(iA),sB(iB))
                error('dimensions to contract should be equal')
            end

            dimsA=setdiff(1:ndims(A),iA);
            dimsB=setdiff(1:ndims(B),iB);


            A=permute(A, [dimsA iA]);
            B=permute(B, [iB dimsB]);

            A=reshape(A, [], prod(sA(iA)));
            B=reshape(B, prod(sB(iB)), []);

            C=A*B;
            output_shape=[sA(dimsA),sB(dimsB)];
            if length(output_shape)>1
                C=reshape(C,[sA(dimsA),sB(dimsB)]);
                if final_permutation
                    C=permute(C,final_permutation);
                end
            end

                function [iA, iB, final_permutation]=parse(s, sA, sB)
                    msg='argument should be a string of the form ''ijk,kjl->il''';
                    if ~ischar(s)
                        error(msg)
                    end

                    %assert that every index appear exactly twice
                    ss=join(split(s,{',','->'}));
                    ss=ss{1};
                    for i=1:length(ss)
                        if length(find(ss==ss(i)))~=2
                            error(['problem with index %s. '...
                                   'Each index should appear exactly twice'], ss(i))
                        end
                    end

                    %split input and output indices
                    s=split(s,'->');
                    if length(s)~=2 
                        error(msg)
                    end

                    %split input indices
                    in=s{1};
                    out=s{2};
                    in=split(in,',');
                    if length(in)~=2
                        error(msg)
                    end
                    inA=in{1};
                    inB=in{2};
                    if length(inA)~=length(sA)
                        error(['''%s'' has %d dimensions while the '...
                            'first argument has %d'],inA, length(inA), length(sA))
                    end
                    if length(inB)~=length(sB)
                        error(['''%s'' has %d dimensions while the '...
                            'second argument has %d'],inB, length(inB), length(sB))
                    end
                    if length(unique(inA))~=length(inA)
                        error('''%s'' has a double index',inA)
                    end
                    if length(unique(inB))~=length(inB)
                        error('''%s'' has a double index',inB)
                    end
                    if length(unique(out))~=length(out)
                        error('''%s'' has a double index',out)
                    end

                    final_permutation=[];
                    iA=[];
                    iB=[];
                    for i=1:length(inA)
                        j=find(inB==inA(i));
                        if isempty(j)   % i is an output index
                            j=find(out==inA(i));
                            final_permutation(end+1)=j;%#ok<AGROW>
                        else            % i is contracted
                            iA(end+1)=i; %#ok<AGROW>
                            iB(end+1)=j; %#ok<AGROW>
                        end
                    end
                    for i=1:length(inB)
                        j=find(inB(i)==out);
                        if ~isempty(j)   % i is an output index
                            final_permutation(end+1)=j;%#ok<AGROW>
                        end
                    end
                    [~, final_permutation]=sort(final_permutation);
                end
        end
    end
    
end