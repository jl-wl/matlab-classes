classdef BlockMBO < handle
    % Many-Body Operator in Blocks according to quantum numbers
    
    properties
        QNumLab % quantum number labels (string array)
        Basis   % basis with quantum numbers and dimensions (struct('qns','dim') array )
        BlkIndr % block row indices (integer array)
        BlkIndc % block column indices (integer array)
        BlkMats % block matrices (cell array)
        QNIndNP % index in QN for fermion number of parity (integer)
    end
    
    methods
        function op = BlockMBO(qnlab, bas, ir, ic, bm, preset)
            arguments
                qnlab   string = []
                bas     struct = struct([])
                ir      double = []
                ic      double = []
                bm      cell = {}
                preset.name {mustBeMember(preset.name, ...
                            {'cu','cd',...
                             'sx','sy','sz','s0'})}
            end
            
            if isfield(preset,'name')
                switch preset.name
                    case 'cu' % c_{up}
                        op.QNumLab = ["N","Sz"];
                        op.Basis = struct('qns',[0,0; 1,1],'dim',[1; 1]);
                        op.BlkIndr = 1;
                        op.BlkIndc = 2;
                        op.BlkMats = {1};
                        op.QNIndNP = 1;
                    case 'cd' % c_{down}
                        op.QNumLab = ["N","Sz"];
                        op.Basis = struct('qns',[0,0; 1,-1],'dim',[1; 1]);
                        op.BlkIndr = 1;
                        op.BlkIndc = 2;
                        op.BlkMats = {1};
                        op.QNIndNP = 1;
                    case 'sx'
                        op.QNumLab = ["N","Sz"];
                        op.Basis = struct('qns',[1,-1; 1,1],'dim',[1; 1]);
                        op.BlkIndr = [1; 2];
                        op.BlkIndc = [2; 1];
                        op.BlkMats = {1; 1};
                        op.QNIndNP = [];
                    case 'sy'
                        op.QNumLab = ["N","Sz"];
                        op.Basis = struct('qns',[1,-1; 1,1],'dim',[1; 1]);
                        op.BlkIndr = [1; 2];
                        op.BlkIndc = [2; 1];
                        op.BlkMats = {1i; -1i};
                        op.QNIndNP = [];
                    case 'sz'
                        op.QNumLab = ["N","Sz"];
                        op.Basis = struct('qns',[1,-1; 1,1],'dim',[1; 1]);
                        op.BlkIndr = [1; 2];
                        op.BlkIndc = [1; 2];
                        op.BlkMats = {-1; 1};
                        op.QNIndNP = [];
                    case 's0'
                        op.QNumLab = ["N","Sz"];
                        op.Basis = struct('qns',[1,-1; 1,1],'dim',[1; 1]);
                        op.BlkIndr = [1; 2];
                        op.BlkIndc = [1; 2];
                        op.BlkMats = {1; 1};
                        op.QNIndNP = [];
                end
                return
            end
            
            if ~isempty(qnlab)
                arrayfun(@(b)assert(size(b.qns,2)==length(qnlab)),bas)
            end
            op.QNumLab = qnlab;
            if ~isfield(bas,'dim') % default dimension = 1 for each QN
                for ii = 1:length(bas)
                    bas(ii).dim = ones(size(bas(ii).qns,1),1);
                end
            end
            op.Basis = bas;
            
            if ~isempty(ir)
                assert(length(bas)==size(ir,2))
            end
            op.BlkIndr = ir;
            
            if isempty(ic)
                op.BlkIndc = ir;
            else
                assert(length(bas)==size(ic,2))
                assert(size(ic,1)==size(ir,1))
                op.BlkIndc = ic;
            end
            
            if ~isempty(bm)
                len = length(bm);
                assert(len==size(ir,1))
                % the following lines are commented out for efficiency
%                 dimr = cell2mat(arrayfun(@(ib)bas(ib).dim(ir(:,ib)),1:length(bas),'UniformOutput',false));
%                 dimr = prod(dimr,2);
%                 dimc = cell2mat(arrayfun(@(ib)bas(ib).dim(ic(:,ib)),1:length(bas),'UniformOutput',false));
%                 dimc = prod(dimc,2);
%                 arrayfun(@(ii)assert(isequal(size(bm{ii}),[dimr(ii),dimc(ii)])),1:len)
            end
            op.BlkMats = reshape(bm,[],1);
            
            op.QNIndNP = find(matches(op.QNumLab,["N","P"],'IgnoreCase',true),1);
        end
        
        function op = addBlock(op, ir, ic, bm) % caution: NOT efficient
            lb = length(op.Basis);
            assert(lb==size(ir,2))
            op.BlkIndr = [op.BlkIndr; ir];
            if isempty(ic)
                op.BlkIndc = [op.BlkIndc; ir];
            else
                assert(lb==size(ic,2))
                op.BlkIndc = [op.BlkIndc; ic];
            end
            dimr = prod(arrayfun(@(ib)op.Basis(ib).dim(ir(ib)),1:lb));
            dimc = prod(arrayfun(@(ib)op.Basis(ib).dim(ic(ib)),1:lb));
            assert(isequal(size(bm),[dimr,dimc]))
            op.BlkMats = [op.BlkMats; {bm}];
        end
        
        function op = select(op, lst)
            op.BlkIndr = op.BlkIndr(lst,:);
            op.BlkIndc = op.BlkIndc(lst,:);
            op.BlkMats = op.BlkMats(lst);
        end
        
        function l = len(op)
            l = size(op.BlkIndr,1);
        end
        
        function qn = ind2qn(op, ind, take_sum)
            qn = zeros(length(ind),length(op.QNumLab));
            for ii = 1:length(ind)
                qn(ii,:) = op.Basis(ii).qns(ind(ii),:);
            end
            if nargin>2 && take_sum
                qn = sum(qn,1);
            end
        end
        
        function qn = qnr(op, ib, take_sum)
            if nargin<3
                take_sum = false;
            end
            qn = op.ind2qn(op.BlkIndr(ib,:),take_sum);
        end
        
        function qn = qnc(op, ib, take_sum)
            if nargin<3
                take_sum = false;
            end
            qn = op.ind2qn(op.BlkIndc(ib,:),take_sum);
        end
        
        function tf = isFermion(op)
            tf = false;
            indnp = op.QNIndNP; % index for fermion number of parity
            if ~isempty(indnp)
                d = op.qnr(1,true)-op.qnc(1,true);
                if mod(d(indnp),2)
                    tf = true;
                end
            end
        end

        function checkDims(op)
        % to check dimension consistency
            [inda, dima] = BlockMBO.allBlocks(op.Basis);
            fdim = @(ii) dima([Utility.findvec(op.BlkIndr(ii,:),inda),Utility.findvec(op.BlkIndc(ii,:),inda)]).';
            arrayfun(@(ii)assert(isequal(size(op.BlkMats{ii}),fdim(ii))),1:op.len)
        end
        
        function tf = hasSameBas(o1, o2)
            if isequal(o1.QNumLab,o2.QNumLab) && isequal(o1.Basis,o2.Basis)
                tf = true;
            else
                tf = false;
            end
        end
        
        function tf = isDiagonal(op, drop0)
        % drop0 is a flag to drop all-zero blocks:
        % if drop0==1, drop only off-diagonal ones
        % if drop0==2, drop both diagonal and off-diagonal ones
            op.gather; % op is in gathered form after testing diagonality
            len = op.len;
            keeplist = true(len,1);
            tf = true;
            for ii = 1:len
                if ~isequal(op.BlkIndr(ii,:),op.BlkIndc(ii,:))
                    if norm(op.BlkMats{ii},'fro')>1e-15
                        tf = false;
                        return
                    else
                        if nargin>1 && drop0
                            keeplist(ii) = false;
                        end
                    end
                else
                    if norm(op.BlkMats{ii},'fro')<1e-15 && nargin>1 && drop0>1
                        keeplist(ii) = false;
                    end
                end
            end
            op.select(keeplist);
        end
        
        function tf = isZero(op)
            tf = false;
            if op.isDiagonal(2)
                if isempty(op.BlkMats)
                    tf = true;
                end
            end
        end
        
        function tf = isEye(op)
            tf = false;
            d = op-op.id;
            if d.isZero
                tf = true;
            end
        end
        
        function op = gather(op)
        % to gather blocks with same QNums (with sorting)
            [blkinds, ~, sortinds] = unique([op.BlkIndr,op.BlkIndc],'rows');
            
            len = size(blkinds,1);
            bm = cell(len,1);
            for ii = 1:len
                bmii = op.BlkMats(sortinds==ii);
                if length(bmii)>1
                    bm{ii} = reshape(sum(reshape(cell2mat(reshape(bmii,1,[])),[],length(bmii)),2),size(bmii{1}));
                else
                    bm(ii) = bmii;
                end
            end
            op.BlkIndr = blkinds(:,1:end/2);
            op.BlkIndc = blkinds(:,1+end/2:end);
            op.BlkMats = bm;
        end
        
        function [opn, rebas] = pack(op, rebas)
        % to pack op (as tensor product) wrt combined QNums with (optional) basis
        % rebas : basis transformation for packed op. It is a struct with fields
        %   - qnn : (opn.Basis={qnsn}) new QNumbers (numeric array)
        %   - qn0 : (qns0=op.Basis) old QNumber basis (cell array)
        %   - ind : indices in old basis for each new QN (cell array)
        %   - dim : dimensions of original blocks corresponding to ind0 (cell array)
        %   - tfm : (optional) additional transformation matrices (cell array)
            if ~isempty(op.QNIndNP) && strcmpi(op.QNumLab(op.QNIndNP),"P")
                indp = op.QNIndNP;
            else
                indp = [];
            end
            
            % new basis
            if nargin<2
                rebas = BlockMBO.packQN(op.Basis, indp);
            end
            
            % block matrices in new basis
            l0 = op.len;
            ir = zeros(l0,1);
            ic = zeros(l0,1);
            bmn = cell(l0,1);
            for ii = 1:l0
                qnr = op.ind2qn(op.BlkIndr(ii,:),true);
                qnc = op.ind2qn(op.BlkIndc(ii,:),true);
                if ~isempty(indp)
                    qnr(indp) = mod(qnr(indp),2);
                    qnc(indp) = mod(qnc(indp),2);
                end
                irii = Utility.findvec(qnr,rebas.qnn);
                icii = Utility.findvec(qnc,rebas.qnn);
                ir(ii) = irii;
                ic(ii) = icii;
                
                cslr = [0; reshape(cumsum(rebas.dim{irii}),[],1)];
                indr = Utility.findvec(op.BlkIndr(ii,:),rebas.ind{irii});
                cslc = [0; reshape(cumsum(rebas.dim{icii}),[],1)];
                indc = Utility.findvec(op.BlkIndc(ii,:),rebas.ind{icii});
                
                bmnii = sparse(cslr(end),cslc(end));
                bmnii(cslr(indr)+1:cslr(indr+1),cslc(indc)+1:cslc(indc+1)) = op.BlkMats{ii}; %#ok<SPRIX>
                if isfield(rebas,'tfm')
                    tfmr = rebas.tfm{irii};
                    tfmc = rebas.tfm{icii};
                    if isempty(tfmr) || isempty(tfmc)
                    % treat empty transformation as ignored HS
                        continue
                    end
                    bmn{ii} = tfmr'*bmnii*tfmc;
                else
                    bmn{ii} = bmnii;
                end
            end
            
            keepblk = ~cellfun(@isempty,bmn);
            [basn, irn, icn] = BlockMBO.transformBas(rebas, ir(keepblk), ic(keepblk));
            opn = BlockMBO(op.QNumLab, basn, irn, icn, bmn(keepblk)).gather;
        end
        
    % --- computational routines ---

        function [ens, rebas, opd, varargout] = diag(op, cutoff, varargin)
        % to diagonalize block-diagonal op with cutoff, and transform other
        % ops in varargin to the diagonalized basis
            if ~isempty(varargin)
                cellfun(@(c)assert(isa(c,'BlockMBO')),varargin)
                cellfun(@(c)assert(op.hasSameBas(c)),varargin)
            end
            
            [opd, rebas] = op.pack;
            if ~opd.isDiagonal(2)
                error('operator to be diagonalized does not preserve symmetry!')
            end
            
            nb = size(opd.Basis.qns,1);
            ens = cell(nb,1);
            tfms = cell(nb,1);
            for ii = 1:nb
                ind = find(opd.BlkIndr==ii);
                if isempty(ind)
                    ens{ii} = zeros(opd.Basis.dim(ii),1);
                	tfms{ii} = speye(opd.Basis.dim(ii));
                else
                    [V, eii] = Utility.eigen(full(opd.BlkMats{ind}));
                    opd.BlkMats{ind} = Utility.spd(eii);
                    ens{ii} = eii;
                	tfms{ii} = V;
                end
            end
            rebas.tfm = tfms;
            
            if nargin>1 && ~isempty(cutoff)
                keepnum = Utility.qfield(cutoff,'number',100); % number of eigenstates to keep
                tolerance = Utility.qfield(cutoff,'tolerance',1e-12); % tolerance for degeneracy
                offset = Utility.qfield(cutoff,'offset',true); % offset smallest eigenvalue to zero
                rescale = Utility.qfield(cutoff,'rescale',1); % rescale eigenvalues
                outfull = Utility.qfield(cutoff,'outfull',false); % output {cut, uncut} of ens, rebas

                if outfull
                    ensf = ens;
                    rebasf = rebas;
                end
    
                if keepnum<sum(cellfun(@length,ens))
                    [ensall, sortinds] = sort(cell2mat(ens));
                    ensall = ensall(ensall<ensall(cutoff.number)+tolerance);
                    sortinds = sortinds(1:length(ensall));
        
                    imax = cumsum([0;cellfun(@length,ens)]);
                    bms = cell(nb,1);
                    for ii = 1:nb
                        ncutii = nnz((sortinds>imax(ii))&(sortinds<=imax(ii+1)));
                        if ncutii>0
                            ens{ii} = ens{ii}(1:ncutii);
                            rebas.tfm{ii} = rebas.tfm{ii}(:,1:ncutii);
                            if offset
                                bms{ii} = Utility.spd(rescale*(ens{ii}-ensall(1)));
                            else
                                bms{ii} = Utility.spd(rescale*ens{ii});
                            end
                        else % discard
                            ens{ii} = [];
                            rebas.tfm{ii} = [];
                            bms{ii} = [];
                        end
                    end
                    opd.BlkIndr = (1:nb)';
                    opd.BlkIndc = (1:nb)';
                    opd.BlkMats = bms;
                    keepblk = ~cellfun(@isempty,bms);
                    opd.select(keepblk);
                    if ~outfull
                        ens = ens(keepblk);
                    end
                end
            else
                outfull = false;
            end

            [opd.Basis, opd.BlkIndr, opd.BlkIndc] = ...
                BlockMBO.transformBas(rebas, opd.BlkIndr, opd.BlkIndc);
            varargout = cellfun(@(o)o.pack(rebas),varargin,'UniformOutput',false);

            if outfull
                ens = {ens, ensf};
                rebas = {rebas, rebasf};
            end
        end

        function c = corr(ham, op1, op2, beta, om)
        % to compute 2-point correlation functions [greater, lesser]
        % lesser has no bosonic/fermionic sign
            [ens, ~, ~, op1, op2] = ham.diag([], op1, op2);
            c = BlockMBO.correlation(op1, op2, ens, beta, om);
        end
        
        function opeye = id(op) % identity operator
            [inda, dima] = BlockMBO.allBlocks(op.Basis);
            opeye = BlockMBO(op.QNumLab, op.Basis, inda, inda,...
                        arrayfun(@(d)speye(d),dima,'UniformOutput',false));
        end
        
        function opP = P(op) % fermion parity operator
            indnp = op.QNIndNP; % index for fermion number of parity
            if isempty(indnp)
                error('parity cannot be determined')
            end
            opP = op.id;
            qnall = arrayfun(@(ii)opP.qnr(ii,true),(1:opP.len)','UniformOutput',false);
            opP.BlkMats = cellfun(@(qn,a)(-1)^qn(indnp)*a,qnall,opP.BlkMats,'UniformOutput',false);
        end
        
        function ope = extend(op, opls, oprs)
        % to extend HS of op by including HS of ops on the
        % left(~outer HS)/right(~inner HS)
        % opls and oprs are both cell array of BlockMBOs
            assert(iscell(opls) && iscell(oprs))
            ope = op;
            if ~ope.isZero && ope.isFermion
                for ii = length(opls):-1:1
                    if opls{ii}.isFermion
                        ope = opls{ii}.P^ope;
                    else
                        ope = opls{ii}.id^ope;
                    end
                end
            else
                for ii = length(opls):-1:1
                    ope = opls{ii}.id^ope;
                end
            end
            for ii = 1:length(oprs)
                ope = ope^oprs{ii}.id;
            end
        end

        function opn = blkwise(op, func, diag_only)
        % to apply func to each block (after gathering)
        % if diag_only (default: false), then require op to be block-diagonal
        % and apply func to both zero and non-zero diagonal blocks;
        % otherwise apply func to each stored (diagonal or not) block in op
            if nargin<2 || isempty(func)
                opn = op;
                return
            end
            if nargin>2 && diag_only
                if ~op.isDiagonal(1)
                    error('operator must be block-diagonal!')
                end
                [inda, dima] = BlockMBO.allBlocks(op.Basis);
                nb = size(inda,1);
                bms = cell(nb,1);
                for ii = 1:nb
                    bi = Utility.findvec(inda(ii,:),op.BlkIndr);
                    if isempty(bi)
                    	bm0 = zeros(dima(ii));
                    else
                    	bm0 = op.BlkMats{bi};
                    end
                    if nargin(func)==3
                        bms{ii} = func(bm0, inda(ii,:), inda(ii,:));
                    elseif nargin(func)==2
                        bms{ii} = func(bm0, inda(ii,:));
                    else
                        bms{ii} = func(bm0);
                    end
                end
                opn = BlockMBO(op.QNumLab, op.Basis, inda, inda, bms);
            else
                op.gather;
                if nargin(func)==3
                    opn = BlockMBO(op.QNumLab, op.Basis, op.BlkIndr, op.BlkIndc, ...
                        cellfun(func,op.BlkMats,num2cell(op.BlkIndr,2),num2cell(op.BlkIndc,2),'UniformOutput',false));
                else
                    opn = BlockMBO(op.QNumLab, op.Basis, op.BlkIndr, op.BlkIndc, ...
                        cellfun(func,op.BlkMats,'UniformOutput',false));
                end
            end
            opn.checkDims
        end

        function s = total(op)
        % to sum over all matrix elements
            s = sum(cellfun(@(m)full(sum(m,'all')), op.BlkMats));
        end
            
    % --- reloaded operators ---

        function p = times(o1, o2)
            assert(o1.hasSameBas(o2))
            if o1.isZero || o2.isZero
                p = o1.null;
                return
            end
            
            ll = o1.len*o2.len;
            ir = zeros(ll,length(o1.Basis));
            ic = zeros(ll,length(o1.Basis));
            bms = cell(ll,1);
            count = 0;
            blkinds1 = [o1.BlkIndr,o1.BlkIndc];
            blkinds2 = [o2.BlkIndr,o2.BlkIndc];
            for ii = 1:o1.len
                for jj = 1:o2.len
                    if isequal(blkinds1(ii,:),blkinds2(jj,:))
                        count = count+1;
                        ir(count,:) = o1.BlkIndr(ii,:);
                        ic(count,:) = o1.BlkIndc(ii,:);
                        bms{count} = o1.BlkMats{ii}.*o2.BlkMats{jj};
                    end
                end
            end
            p = BlockMBO(o1.QNumLab, o1.Basis, ir(1:count,:), ic(1:count,:), bms(1:count));
            p.gather;
        end

        function p = mtimes(o1, o2)
        % product of operators in same HS, including number*operator
            if isnumeric(o1)
                p = BlockMBO(o2.QNumLab, o2.Basis, o2.BlkIndr, o2.BlkIndc, ...
                    cellfun(@(a)o1*a,o2.BlkMats,'UniformOutput',false));
                return
            elseif isnumeric(o2)
                p = BlockMBO(o1.QNumLab, o1.Basis, o1.BlkIndr, o1.BlkIndc, ...
                    cellfun(@(a)o2*a,o1.BlkMats,'UniformOutput',false));
                return
            end
            
            assert(o1.hasSameBas(o2))
            if o1.isZero || o2.isZero
                p = o1.null;
                return
            end
            
            ll = o1.len*o2.len;
            ir = zeros(ll,length(o1.Basis));
            ic = zeros(ll,length(o1.Basis));
            bms = cell(ll,1);
            count = 0;
            for ii = 1:o1.len
                for jj = 1:o2.len
                    if isequal(o1.BlkIndc(ii,:),o2.BlkIndr(jj,:))
                        count = count+1;
                        ir(count,:) = o1.BlkIndr(ii,:);
                        ic(count,:) = o2.BlkIndc(jj,:);
                        bms{count} = o1.BlkMats{ii}*o2.BlkMats{jj};
                    end
                end
            end
            p = BlockMBO(o1.QNumLab, o1.Basis, ir(1:count,:), ic(1:count,:), bms(1:count));
            p.gather;
        end
            
        function oo = mpower(o2, o1)
        % if o1 or o2 is numeric scalar and the other is block-diagonal op,
        % then perform normal mpower;
        % otherwise used as short-hand (^) for operator tensor product,
        % where o2^o1 ~ kron(o2,o1), thus HS1 is interior to HS2
            if isnumeric(o1)
                if ~o2.isDiagonal(1)
                    error('operator must be block-diagonal!')
                end
                oo = BlockMBO(o2.QNumLab, o2.Basis, o2.BlkIndr, o2.BlkIndc, ...
                    cellfun(@(a)a^o1,o2.BlkMats,'UniformOutput',false));
                return
            elseif isnumeric(o2)
                if ~o1.isDiagonal(1)
                    error('operator must be block-diagonal!')
                end
                oo = o1.blkwise(@(m)o2^m, true);
                return
            end

            % tensor product if both are operators
            assert(isequal(o1.QNumLab,o2.QNumLab))
            bas = [o1.Basis, o2.Basis];
            if o1.isZero || o2.isZero
                oo = BlockMBO(o1.QNumLab, bas);
                return
            end
            
            indnp = o1.QNIndNP; % index for fermion number of parity
            if ~isempty(indnp)
                d1 = o1.qnr(1,true)-o1.qnc(1,true);
                d2 = o2.qnr(1,true)-o2.qnc(1,true);
                fex = mod(d1(indnp)*d2(indnp),2); % 1: fermionic exchange
            else
                fex = 0;
            end
            
            ll = o1.len*o2.len;
            
            ir = zeros(ll,length(bas));
            ic = zeros(ll,length(bas));
            bms = cell(ll,1);
            count = 0;
            for ii = 1:o2.len
                qnc2 = o2.qnc(ii,true);
                bm2 = o2.BlkMats{ii};
                for jj = 1:o1.len
                    bm = kron(bm2, o1.BlkMats{jj});
                    if fex
                        bm = (-1)^qnc2(indnp)*bm;
                    end
                    count = count+1;
                    ir(count,:) = [o1.BlkIndr(jj,:),o2.BlkIndr(ii,:)];
                    ic(count,:) = [o1.BlkIndc(jj,:),o2.BlkIndc(ii,:)];
                    bms{count} = bm;
                end
            end
            oo = BlockMBO(o1.QNumLab, bas, ir, ic, bms);
        end
        
        function opn = expm(op)
        % to do block-diagonal-matrix exponential (TODO: generalize)
            opn = op.blkwise(@(m)expm(m), true);
        end

        function s = plus(o1, o2)
            assert(o1.hasSameBas(o2))
            s = BlockMBO(o1.QNumLab, o1.Basis,...
                         [o1.BlkIndr; o2.BlkIndr],...
                         [o1.BlkIndc; o2.BlkIndc],...
                         [o1.BlkMats; o2.BlkMats]);
            s.gather;
        end
        
        function po = uplus(op)
            po = op;
        end
        
        function mop = uminus(op)
            mop = BlockMBO(op.QNumLab, op.Basis, op.BlkIndr, op.BlkIndc, ...
                cellfun(@(a)-a,op.BlkMats,'UniformOutput',false));
        end
        
        function d = minus(o1, o2)
            d = o1+(-o2);
        end
        
        function opt = transpose(op)
            opt = BlockMBO(op.QNumLab, op.Basis);
            opt.BlkIndr = op.BlkIndc;
            opt.BlkIndc = op.BlkIndr;
            opt.BlkMats = cellfun(@transpose,op.BlkMats,'UniformOutput',false);
        end
        
        function opct = ctranspose(op)
            opct = BlockMBO(op.QNumLab, op.Basis);
            opct.BlkIndr = op.BlkIndc;
            opct.BlkIndc = op.BlkIndr;
            opct.BlkMats = cellfun(@ctranspose,op.BlkMats,'UniformOutput',false);
        end

        function t = trace(op)
            t = 0;
            for ii = 1:op.len
                if isequal(op.BlkIndr(ii,:),op.BlkIndc(ii,:))
                    t = t + trace(op.BlkMats{ii});
                end
            end
        end
        
        function opn = null(op)
            opn = BlockMBO(op.QNumLab, op.Basis);
        end
    end
    
    methods (Static)
        function [inda, dima, qnsa] = allBlocks(bas)
        % indices, dimensions, quantum numbers of all blocks
            qn0 = {bas.qns};
            ls = cellfun(@(a)size(a,1),qn0);
            inda = Utility.ndgridVecs(ones(size(ls)),ls);
            if nargout>1
            dima = prod(arrayfun(@(ind,ii)bas(ii).dim(ind),...
                                 inda,repmat(1:length(ls),prod(ls),1)),...
                        2);
            if nargout>2
            qnsa = reshape(cell2mat(...
                       arrayfun(@(ind,ii)bas(ii).qns(ind,:),...
                                inda,repmat(1:length(ls),prod(ls),1),...
                                'UniformOutput',false)),...
                       prod(ls),[],length(ls));
            end
            end
        end

        function rebas = packQN(bas0, indp)
        % to pack QNumbers in old basis (as op.Basis) to form new basis
        % rebas : basis transformation information. It is a struct with fields
        %   - qnn : (opn.Basis={qnsn}) new QNumbers (numeric array)
        %   - qn0 : (qns0={op.Basis.qns}) old QNumber basis (cell array)
        %   - ind : indices in old basis for each new QN (cell array)
        %   - dim : dimensions of original blocks corresponding to ind0 (cell array)
            if length(bas0)==1
%                 warning('basis does not need to be packed')
                rebas.qnn = bas0.qns;
                rebas.qn0 = {bas0.qns};
                rebas.ind = num2cell((1:size(bas0.qns,1))');
                rebas.dim = num2cell(bas0.dim);
                return
            end
            
            [inda, dima, qnsa] = BlockMBO.allBlocks(bas0);
                    
            [qnn, ~, indn] = unique(sum(qnsa,3),'rows');
            if nargin>1 && ~isempty(indp)
                qnn(:,indp) = mod(qnn(:,indp),2);
            end

            rebas.qnn = qnn;
            rebas.qn0 = {bas0.qns};
            rebas.ind = arrayfun(@(ii)inda(indn==ii,:),(1:size(qnn,1))','UniformOutput',false);
            rebas.dim = arrayfun(@(ii)dima(indn==ii),(1:size(qnn,1))','UniformOutput',false);
        end
        
        function [basn, irn, icn] = transformBas(rebas, ir, ic)
        % to form new basis and corresponding block-indices according to rebas
            if isfield(rebas,'tfm')
                keepqn = find(~cellfun(@isempty,rebas.tfm));
                basn.qns = rebas.qnn(keepqn,:);
                basn.dim = nonzeros(cellfun(@(m)size(m,2),rebas.tfm));
                indqn = zeros(size(rebas.qnn,1),1);
                indqn(keepqn) = (1:length(keepqn))';
                irn = indqn(ir);
                icn = indqn(ic);
            else
                basn.qns = rebas.qnn;
                basn.dim = cellfun(@sum,rebas.dim);
                irn = ir;
                icn = ic;
            end
        end
        
        function wfs = rebas2wfs(rebas)
            assert(isfield(rebas,'tfm'))
            wfs.qns = rebas.qnn;
            wfs.amp = rebas.tfm;
            
            ind = rebas.ind;
            bas = cell(length(ind),1);
            for ii = 1:length(ind)
                indii = ind{ii};
                basii = cell(size(indii));
                for jj = 1:size(indii,1)
                    basii(jj,:) = arrayfun(@(kk)rebas.qn0{kk}(indii(jj,kk),:),...
                                    1:length(rebas.qn0),'UniformOutput',false);
                end
                bas{ii} = basii;
            end
            wfs.bas = bas;
        end

        function c = correlation(op1, op2, ens, beta, om, mask)
        % to compute [greater, lesser] correlation functions
        % op1 and op2 are in the basis of diagonalized Hamiltonian
        % lesser has no bosonic/fermionic sign
            enmin = min(cell2mat(ens));
            ens = cellfun(@(a)a-enmin,ens,'UniformOutput',false);
            mtimes_w = @(m,ir,ic) m.*exp(-beta*ens{ic}.');
            w_mtimes = @(m,ir,ic) exp(-beta*ens{ir}).*m;
            times_g0 = @(m,ir,ic) m./(om+ens{ir}-ens{ic}.');
            o2g = op2.transpose.blkwise(times_g0);
            cgtr = op1.blkwise(w_mtimes).*o2g;
            clsr = op1.blkwise(mtimes_w).*o2g;
            Z = sum(cellfun(@(a)sum(exp(-beta*a)),ens));
            if nargin<6
                c = [cgtr.total, clsr.total]/Z;
            else
                c = [cgtr.blkwise(mask).total, clsr.blkwise(mask).total]/Z;
            end
        end
        
        function sv = svec(type)
            switch type
                case 'c'
                    cu = BlockMBO('name','cu');
                    cd = BlockMBO('name','cd');
                    sv = {-cd^(cu')+cd'^cu, 1i*(cd^(cu')+(cd')^cu),...
                          cd.id^(cu'*cu)-(cd'*cd)^cu.id};
                case 's'
                    sv = {BlockMBO('name','sx'), BlockMBO('name','sy'),...
                          BlockMBO('name','sz')};
            end
        end
        
        function test
            cu = BlockMBO('name','cu');
            cd = BlockMBO('name','cd');
            sx = BlockMBO('name','sx');
            sy = BlockMBO('name','sy');
            sz = BlockMBO('name','sz');
            
            % test Fermion statistics
            cuex = cu.extend({cd},{});
            cdex = cd.extend({},{cu});
            cellfun(@(op)assert(op.isZero),...
                    {cuex*cuex, cdex*cdex,...
                     cuex*cdex+cdex*cuex,...
                     cuex'*cdex+cdex*cuex'});
            cellfun(@(op)assert(op.isEye),...
                    {cuex'*cuex+cuex*cuex',...
                     cdex'*cdex+cdex*cdex'});
                 
            % test commutation relation of Pauli matrices
            cellfun(@(op)assert(op.isZero),...
                    {sx*sy+sy*sx,sy*sz+sz*sy,sz*sx+sx*sz});
            cellfun(@(op)assert(op.isEye),...
                    {sx*sx,sy*sy,sz*sz});
            cellfun(@(op)assert(op.isZero),...
                    {sx*sy-1i*sz,sy*sz-1i*sx,sz*sx-1i*sy});
            
            disp('all tests passed')
        end
    end
    
end