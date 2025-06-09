classdef TBHamiltonian
    % Tight-Binding Hamiltonian
    % TODO: handle class?
    % TODO: setUnitCell (update DoF.Site with labels)
    
    properties
        Vectors % [NVecs, NDims]; R = to - from
        Matrices % sparse: [NDMat*NDMat, NVecs]
        DoF % degrees of freedom (struct array)
        UnitCell
        Partition % indices for each site (cell array of column vectors)
        NDims % spacial dimension
        NVecs % number of vectors in Rs
        NDMat % number of dimensions of HR
    end
    
    methods
        function h = TBHamiltonian(Rs, HR, dof, uc)
        % when both UnitCell and DoF are nonempty, Partition is always defined
            [h.Vectors, h.Matrices] = TBHamiltonian.combineSameR(Rs,HR);
            h.NDims = size(h.Vectors,2);
            h.NVecs = size(h.Vectors,1);
            h.NDMat = sqrt(size(h.Matrices,1));
            if nargin>2 && ~isempty(dof)
                if iscell(dof)
                    [h.DoF, h.Partition] = TBHamiltonian.convertDoF(dof);
                elseif isstruct(dof)
                    h.DoF = dof;
                    h.Partition = [];
                else
                    error('Invalid dof input!')
                end
            else
                h.DoF = [];
                h.Partition = [];
            end
            if ~isempty(h.DoF) && length(h.DoF)~=h.NDMat
                error('Hamiltonian and DoF do not match!')
            end
            if nargin>3 && ~isempty(uc)
                if uc.NDims~=h.NDims
                    error('Incompatible number of dimensions!')
                end
                if iscell(dof)
                    labels = cellfun(@(c)c{1},dof,'UniformOutput',false);
                    if uc.NSites~=length(dof)
                        error('Incompatible number of sites!')
                    elseif isempty(uc.Labels)
                        uc.Labels = labels;
                    elseif ~isempty(setxor(uc.Labels,labels))
                        error('Incompatible labels of sites!')
                    end
                elseif isstruct(dof)
                    if isempty(uc.Labels)
                        h.Partition = {(1:h.NDMat).'};
                    elseif isfield(dof,'Site')
                        labels = unique({dof.Site});
                        if ~isempty(setxor(uc.Labels,labels))
                            error('Incompatible labels of sites!')
                        else
                            h.Partition = cellfun(@(l)TBHamiltonian.findIndices(dof,'Site',l), ...
                                                  uc.Labels,'UniformOutput',false);
                        end
                    elseif isscalar(uc.Labels)
                        s = repelem(uc.Labels,h.NDMat,1);
                        [h.DoF.Site] = s{:};
                        h.Partition = {(1:h.NDMat).'};
                    end
                end
                h.UnitCell = uc;
            else
                h.UnitCell = [];
            end
        end
        
        function dof = getDoF(h, ucrs)
            dof0 = h.DoF;
            if nargin<2
                dof = dof0;
                return
            end
            dH = h.NDMat;
            nuc = size(ucrs,1);
            crs = mat2cell(ucrs,ones(nuc,1)).';
            if isempty(dof0)
                c = repelem(crs,1,dH);
                f = {'UnitCell'};
            else
                c = vertcat(repelem(crs,1,dH),repmat(struct2cell(reshape(dof0,[],1)),1,nuc));
                f = vertcat('UnitCell',fieldnames(dof0));
            end
            dof = cell2struct(c,f);
        end
        
        function [pos, part] = getPosition(h, ucrs, frac_coordinates)
        % positions and partitions for all sites
        % position starts from 0 by default
            if nargin<3
                frac_coordinates = true; % fractional by default!
            end
            uc = h.UnitCell;
            part = h.Partition;
            nd = h.NDims;
            dH = h.NDMat;
            if isempty(uc) || isempty(part)
%                 warning(['No unit cell or partition information, ', ...
%                          'set all positions in one unit cell to origin.'])
                pos = zeros(1,nd);
                part = {(1:dH).'};
            else
                pos = uc.Sites;
            end
            if nargin>1 && ~isempty(ucrs)
                nuc = size(ucrs,1);
                pos = repmat(pos,nuc,1)+repelem(ucrs,size(pos,1),1);
                np = length(part);
                part = cellfun(@(a,b)a+b,repmat(part,nuc,1),...
                           mat2cell(repelem(dH*(0:nuc-1).',np,1),ones(nuc*np,1)),...
                           'UniformOutput',false);
            end
            if ~frac_coordinates
                pos = pos*uc.Basis+repmat(uc.Origin,size(pos,1),1);
            end
        end
        
        function h = addHops(h, Rs, HR)
            Rsn = cat(1,h.Vectors,Rs);
            HRn = cat(2,h.Matrices,sparse(reshape(HR,[],size(Rs,1))));
            [h.Vectors, h.Matrices] = TBHamiltonian.combineSameR(Rsn,HRn);
            h.NVecs = size(h.Vectors,1);
        end
        
        function h = setChemicalPotential(h, mu)
            ind0 = find(~sum(abs(h.Vectors),2));
            if ~isfield(h.DoF,'Nambu')
                a = speye(h.NDMat);
            else
                a = TBHamiltonian.fullOperator(h.DoF,{'Nambu',{'e','h'},[1,0;0,-1]});
            end
            h.Matrices(:,ind0) = h.Matrices(:,ind0)-mu*reshape(a,[],1);
        end
        
        function hn = reconstructUnitCell(h0, new_basis, new_origin, rotation)
            Rs0 = h0.Vectors;
            HR0 = h0.Matrices;
            uc0 = h0.UnitCell;
            dof0 = h0.DoF;
            part0 = h0.Partition;
            nR0 = h0.NVecs;
            dH0 = h0.NDMat;
            
            if nargin<2
                error('Must at least provide new basis!')
            elseif nargin<3
                ucn = uc0.reconstruct(new_basis);
            else
                ucn = uc0.reconstruct(new_basis,new_origin);
            end
            if nargin<4
                rotation = 1;
            end
            
            bas0 = uc0.Basis;
            pos0 = uc0.Sites;
            n0 = uc0.NSites;
            basn = ucn.Basis;
            posn = ucn.Sites;
            n = ucn.NSites;
            m = n/n0;
            
            dHn = m*dH0;
            lp0 = cellfun(@length,part0);
            partn = mat2cell((1:dHn).',repmat(lp0,m,1));
            dofn = repmat(dof0(cell2mat(part0)),m,1);
            labn = repelem(ucn.Labels,cellfun(@length,partn),1);
            [dofn.Site] = labn{:};
            
            posn2 = repelem(posn,n0*nR0,1)...
                    +(repmat((repelem(pos0,nR0,1)+repmat(Rs0,n0,1)),n,1)...
                      -repmat(repelem(pos0,n0*nR0,1),m,1))*bas0/basn;
            Rsn = floor(posn2);
            
            indn = cellfun(@(p)find(all(abs(posn-repmat(p,n,1))<1e-12,2)),...
                           mat2cell(mod(posn2,1),ones(size(posn2,1),1)));
            HRn = arrayfun(@copy_hr,indn.',...
                           repelem(1:n,1,n0*nR0),...
                           repmat(repelem(1:n0,1,nR0),1,n),...
                           repmat(repelem(1:n0,1,n0*nR0),1,m),...
                           repmat(repmat(1:nR0,1,n0),1,n),...
                           'UniformOutput',false);
            HRn = cell2mat(HRn);
            
            ucn = ucn.rotateCoordinates(rotation);
            % rotate only after constructing new Hamiltonian
            
            hn = TBHamiltonian(Rsn,HRn,dofn,ucn);
%             hn.Partition = partn; % to ensure same order after copying
            
            function h = copy_hr(rn, cn, r0, c0, iR0)
                h = sparse(dHn*dHn,1);
                h(sub2ind([dHn,dHn],repmat(partn{rn},lp0(c0),1),...
                                    repelem(partn{cn},lp0(r0),1)))...
                = HR0(sub2ind([dH0,dH0],repmat(part0{r0},lp0(c0),1),...
                                        repelem(part0{c0},lp0(r0),1)),...
                      iR0);
            end
        end

        function hh = hc(h)
            hh = h;
            nR = size(hh.Vectors,1);
            if nR>0
                hh.Vectors = -hh.Vectors;
                d = sqrt(size(hh.Matrices,1));
                for iR = 1:nR
                    m = reshape(hh.Matrices(:,iR),d,d);
                    hh.Matrices(:,iR) = reshape(m',[],1);
                end
            end
        end
 
        function [Rsn, HRn] = transform(h0, trans)
        % transform Hamiltonian according function trans
            Rs0 = h0.Vectors;
            HR0 = h0.Matrices;
            dH0 = h0.NDMat;
            
            cRs0 = num2cell(Rs0,2);
            cHR0 = num2cell(HR0,1);
            cHR0 = cellfun(@(a)reshape(a,dH0,dH0),cHR0,'UniformOutput',false);
            if nargin(trans)==1
                cHRn = cellfun(@(a)trans(a),...
                               cHR0,'UniformOutput',false);
                cRsn = cRs0;
            elseif nargin(trans)==2
                cHRn = cellfun(@(a)trans(a,h0),...
                               cHR0,'UniformOutput',false);
                cRsn = cRs0;
            else
                [cHRn, cRsn] = cellfun(@(a,r)trans(a,r,h0),...
                               cHR0,cRs0.','UniformOutput',false);
            end
            Rsn = cell2mat(reshape(cRsn,[],1));
            HRn = cell2mat(reshape(cHRn,1,[]));
            [Rsn, HRn] = TBHamiltonian.combineSameR(Rsn, HRn);
        end
        
        function varargout = isSymmetry(h0, trans, tol)
        % check if trans is a symmetry for Hamiltonian h0 with tolerance tol
            Rs0 = h0.Vectors;
            HR0 = h0.Matrices;
            [Rsn, HRn] = h0.transform(trans);
            tfR = isequal(Rs0,Rsn);
            if nargin<3
                tol = 1e-12;
            end
            if isequal(size(HR0),size(HRn))
                tfH = ~full(any(Utility.chop(HR0-HRn,tol),'all'));
            else
                tfH = false;
            end
            if nargout==0
                disp(tfR & tfH)
            elseif nargout==1
                varargout{1} = tfR & tfH;
            elseif nargout==2
                varargout{1} = tfR;
                varargout{2} = tfH;
            end
        end
        
        function hn = double(h0, dofn, transform, add_terms)
        % double DoF of Ham, e.g. to add SOC or pairing
        % dofn can be full struct array for the new Hamiltonian or cell
        %   array containing only the new dof, e.g. {'Spin',{'u','d'}} or
        %   {'Nambu',{'e','h'}}
        % transform(H[[,R],DoF]) is a function handle
        % add_terms has the form {R1,mat1,R2,mat2,...}
            Rs0 = h0.Vectors;
            HR0 = h0.Matrices;
            dof0 = h0.DoF;
            uc0 = h0.UnitCell;
            nR0 = h0.NVecs;
            dH0 = h0.NDMat;
            
            cHR0 = mat2cell(HR0,size(HR0,1),ones(nR0,1));
            cHR0 = cellfun(@(a)reshape(a,dH0,dH0),cHR0,'UniformOutput',false);
            if nargin(transform)==1
                cHRn = cellfun(@(a)blkdiag(a,transform(a)),...
                               cHR0,'UniformOutput',false);
            elseif nargin(transform)==2
                cHRn = cellfun(@(a)blkdiag(a,transform(a,dof0)),...
                               cHR0,'UniformOutput',false);
            else
                cHRn = cellfun(@(a,r)blkdiag(a,transform(a,r,dof0)),...
                               cHR0,mat2cell(Rs0,ones(nR0,1)),'UniformOutput',false);
            end
            if nargin>3 && ~isempty(add_terms)
                if ~iscell(add_terms)
                    error('Invalid additional terms!')
                end
                for iR = 1:length(add_terms)/2
                    R = add_terms{2*iR-1};
                    ind = find(~sum(abs(Rs0-repmat(R,nR0,1)),2));
                    if isempty(ind)
                        Rs0 = [Rs0; R];
                        cHRn{end+1} = add_terms{2*iR};
                    else
                        cHRn{ind} = cHRn{ind}+add_terms{2*iR};
                    end
                    mR = find(~cellfun(@(a)norm(a+R),add_terms(1:2:end))); %#ok<EFIND>
                    if isempty(mR) % no opposite R; add conjugate terms
                        ind = find(~sum(abs(Rs0+repmat(R,nR0,1)),2));
                        if isempty(ind)
                            Rs0 = [Rs0; -R];
                            cHRn{end+1} = add_terms{2*iR}';
                        else
                            cHRn{ind} = cHRn{ind}+add_terms{2*iR}';
                        end
                    end
                    nR0 = size(Rs0,1);
                end
            end
            HRn = reshape(cell2mat(cHRn),[],nR0);
            
            if iscell(dofn)
                name = dofn{1};
                vals = dofn{2};
                dofn = repmat(dof0,2,1);
                for ii = 1:dH0
                    dofn(ii).(name) = vals{1};
                    dofn(ii+dH0).(name) = vals{2};
                end
            end
            
            hn = TBHamiltonian(Rs0,HRn,dofn,uc0);
        end
        
        function es = bulkSpectrum(h, ks) % ks are in units of 2pi/a
            Rs = h.Vectors;
            HR = h.Matrices;
            nd = h.NDims;
            dH = h.NDMat;

            if size(ks,2)~=nd
                if size(ks,1)==nd
                    ks = ks.';
                else
                    error('Space dimensions do not match!')
                end
            end

            eigk = @(k) Utility.eigen(full(reshape(HR*exp(-2i*pi*Rs*k'),dH,dH)));

            es = cellfun(eigk,mat2cell(ks,ones(size(ks,1),1)),'UniformOutput',false);
            es = reshape(cell2mat(es),dH,[]);
        end
        
        function [ham, pos, part, dof] = onLattice(h, latps)
        % latps can contain the following parameters:
        %   ns: number of unit cells in each dimension, 0 for infinity
        % 	ks: momentum in units of 2pi/a 
        %   pbcs: periodical boundary conditions (multiplied to hopping terms; default 0)
        %   Bz: magnetic field in units of Tesla
        %   frac: flag for using fractional coordinates in pos
        %   trim: conditions to trim the lattice
        %   ost: onsite terms
        %   rcp: rescaling factor/function for coupling between sites
        % Elements of ks are valid only in infinite dims;
        % When Bz is nonempty, Landau gauge A=(0,Bz*x,0...) is used,
        % the 1st dim has to be finite and its basis vector has to be
        % orthogonal to all other basis vectors.
            if isfield(latps,'ns')
                s = size(latps.ns);
            elseif isfield(latps,'ks')
                s = size(latps.ks);
            else
                error('Must provide lattice size or momentum vector!')
            end

            qfield = @Utility.qfield;
            ns = qfield(latps,'ns',zeros(s));
            ks = qfield(latps,'ks',zeros(s));
            pbcs = qfield(latps,'pbcs',zeros(s));
            Bz = qfield(latps,'Bz',[]);

            nd = h.NDims;
            if length(ns)~=nd || length(ks)~=nd || length(pbcs)~=nd
                error('Space dimensions do not match!')
            end

            if ~isempty(Bz)
                if nd<2
                    error('Magnetic field is valid only for higher than 1-dimension!')
                end
                if ns(1)==0
                    error('The first dimension must be finite to apply Landau gauge!')
                end
            end

            ks = ks.*(ns==0); % momentum only valid for infinite dimensions
            pbcs(ns==0) = 1; % set infinite dimensions to be periodical
            ns(ns==0) = 1; % treat infinite dimensions as single unit cell
            ucrs = Utility.ndgridVecs(zeros(size(ns)),ns-1); % Rs for all unit cells

            hh = h;
            while nd>1
                Rs = hh.Vectors;
                nR = hh.NVecs;
                % separate Rs according to the second to the last dims
                step = find(sum(abs(diff(Rs(:,2:end),1,1)),2)>0.01);
                cinds = mat2cell(1:nR,1,diff([0;step;nR]));
                HR = cell2mat(cellfun( ...
                         @(inds)reshape(hh.extract(inds).ham1D(ns,ks,pbcs,Bz),[],1), ...
                         cinds, 'UniformOutput', false));
                % reduce dimensions by 1
                hh = TBHamiltonian(Rs([step;nR],2:end),HR);
                ns = ns(2:end);
                ks = ks(2:end);
                pbcs = pbcs(2:end);
                if ~isempty(Bz)
                    Bz = [];
                end
                nd = nd-1;
            end
            ham = hh.ham1D(ns,ks,pbcs);
            
            frac = qfield(latps,'frac',true); % fractional by default!
            [pos, part] = h.getPosition(ucrs,frac);
            dof = h.getDoF(ucrs);
            
            % trim lattice
            if isfield(latps,'trim')
                trim = latps.trim;
                % trim is rows of vectors, e.g. in 3D,
                % [a,b,c,d] such that the trim condition is ax+by+cz+d>0
                % each row gives one such trim condition
                % update: trim can be cell array of the above form
                % fulfilling "or" function
                % update: trim can be function handle
                if isa(trim,'function_handle')
                    keep = trim(pos);
                else
                    fkeep = @(t) find(all([pos,ones(size(pos,1),1)]*t.'>0,2));
                    if iscell(trim)
                        ckeep = cellfun(fkeep,trim,'UniformOutput',false);
                        keep = any(cell2mat(reshape(ckeep,1,[])),2);
                    else
                        keep = fkeep(trim);
                    end
                end
                pos = pos(keep,:);
                part = part(keep);
                ind = cat(1,part{:});
                ham = ham(ind,ind); % trimming will reorganize Ham basis
                part = mat2cell((1:length(ind)).',cellfun(@length,part),1);
                dof = dof(ind);
            end
            
            % add onsite terms
            if isfield(latps,'ost')
                ost = latps.ost;
                if ~iscell(ost)
                    ost = {ost};
                end
                for ii = 1:length(ost)
                    ham = TBHamiltonian.addOnsiteTerm(ost{ii}, ham, pos, part, dof);
                    % note here pos is controlled by latps.frac
                end
            end
            
            % rescale coupling between sites
            if isfield(latps,'rcp')
                rcp = latps.rcp;
                if size(rcp,1)==3 && size(rcp,2)==1
                    rcp = rcp';
                end
                for ii = 1:size(rcp,1)
                    r = rcp{ii,1};
                    mask1 = rcp{ii,2};
                    mask2 = rcp{ii,3};
                    ham = TBHamiltonian.rescaleSiteCoupling(r, ham, part, mask1, mask2);
                end
            end
        end
        
        function [spf, ucrs] = surfaceSpectralFun(h, ia, k, ens, eni)
        % spf contains all dof in the surface principle layer,
        % ucrs is the Rs for the unit cells in the surface principle layer.
        % The surface is defined such that: if ia>0 (ia<0), the system
        % is semi-infinite with R(|ia|)>=0 (R(|ia|)<=0).
            if nargin<5
                eni = 1e-5;
            end
            enpts = length(ens);
            [h0, v, maxR] = principalLayer(h,ia,k);
            d = size(h0,1);
            spf = zeros(d,enpts);
            for ien = 1:enpts
                en = ens(ien)+1i*eni;
                tm = ComputePhys.tMat(en*speye(d)-h0,v);
                sgf = inv(en*speye(d)-h0-v'*tm);
                spf(:,ien) = -imag(diag(sgf));
            end
            ucrs = zeros(maxR,h.NDims);
            if ia>0
                ucrs(:,ia) = (0:maxR-1).';
            else
                ucrs(:,-ia) = (1-maxR:0).';
            end
        end
        
        function [sgf, h0, v, maxRa, pos, part, dof] = surfaceGreenFun(h, latps, ens, eni)
        % sgf contains all dof in the surface principle layer.
        % For values in latps.ns:
        %   0 means infinite,
        %   -1 means semi-infinite in positive direction (s>0),
        %   -2 means semi-infinite in negative direction (s<0),
        %   there is one and only one negative value (ns(ia)) in latps.ns.
        % The surface is defined such that:
        %   if s>0, the system is semi-infinite with R(ia)>=0, and
        %   v is from i-th principal layer to (i+1)-th;
        %   if s<0, the system is semi-infinite with R(ia)<=0, and
        %   v is from i-th principal layer to (i-1)-th.
        % The self-energy will always be v'*sgf*v
            if nargin<4
                eni = 1e-5;
            end
            enpts = length(ens);
            [h0, v, maxRa, pos, part, dof] = principalLayer(h,latps);
            d = size(h0,1);
            sgf = zeros(d,d,enpts);
            for ien = 1:enpts
                en = ens(ien)+1i*eni;
                tm = ComputePhys.tMat(en*speye(d)-h0,v);
                sgf(:,:,ien) = inv(en*speye(d)-h0-v'*tm);
            end
        end
        
        function gf = bulkGreenFun(h, latps, ens, eni, kpts)
        % Bulk GF with possible finer k-sampling (inserted through kpts)
        % in order to account for large embedding lattice.
        % Periodical boundary condition by default.
            if isfield(latps,'ns')
                s = size(latps.ns);
            elseif isfield(latps,'ks')
                s = size(latps.ks);
            else
                error('Must provide lattice size or momentum vector!')
            end

            if nargin<4
                eni = 1e-5;
            end
            enpts = length(ens);
            
            if ~isfield(latps,'pbcs')
%                 warning('Convert to periodical boundary condition!')
                latps.pbcs = ones(s);
            end
            ham = h.onLattice(latps);
            d = size(ham,1);
            gf = zeros(d,d,enpts);
            for ien = 1:enpts
                en = ens(ien)+1i*eni;
                gf(:,:,ien) = inv(en*eye(d)-ham);
            end
            
            if nargin<5 || ~isfield(latps,'ns') || ~any(latps.ns>0) % no finer k-sampling
                return
            elseif nnz(latps.pbcs~=0 & latps.ns>0) ~= nnz(latps.ns>0)
                warning('Incompatible open boundary condition with k-sampling insertion!')
            end

            if isscalar(kpts)
                kpts = kpts*ones(s);
            end
            kpts(latps.ns<=0) = 1; % no insertion if not finite size
            dks = Utility.ndgridVecs(zeros(s),kpts-1,1)*diag(1./kpts);

            for idk = 2:size(dks,1) % skip origin
                ham = h.onLattice(setfield(latps,'pbcs',exp(2i*pi*dks(idk,:)))); %#ok<SFLD> 
                for ien = 1:enpts
                    en = ens(ien)+1i*eni;
                    gf(:,:,ien) = gf(:,:,ien) + inv(en*eye(d)-ham);
                end
            end
            gf = gf/size(dks,1);
        end

        function [h0, v, maxRa, pos, part, dof] = principalLayer(h, varargin)
        % Hamiltonian of intra (h0) and inter (v) principal layers.
        % Accept two types of inputs: (h, ia, k) or (h, latps)
        % If (as in old version) input is (h, ia, k), then
        %   ia (s=sign(ia)) is the axis perp to which layers are defined,
        %   k is the momentum (in units of 2pi/a) parallel to the layers.
        %   Note that dim(k) = dim(R).
        % If input is (h, latps), then for values in latps.ns,
        %   0 means infinite,
        %   -1 means semi-infinite in positive direction (s=+1),
        %   -2 means semi-infinite in negative direction (s=-1),
        %   there is one and only one negative value in latps.ns.
        % In output, v is from i-th principal layer to (i+s)-th
            Rs = h.Vectors;
            nd = h.NDims;
            nR = h.NVecs;
            
            if nargin==2
                latps = varargin{1};
                ns = latps.ns;
                if nnz(ns==-1 | ns==-2)~=1
                    error('Invalid setting in latps.ns!')
                end
                ia = find(ns<0,1);
                switch ns(ia)
                    case -1
                        s = 1; % direction of semi-infinity
                    case -2
                        s = -1;
                    otherwise
                        error('Invalid setting in latps.ns!')
                end
                if ~isfield(latps,'ks')
                    latps.ks = zeros(1,nd);
                end
            elseif nargin==3
                latps.ns = zeros(1,nd); % all dimensions are infinite except for ia 
                latps.ks = varargin{2};
                ia = varargin{1};
                s = sign(ia); % direction of semi-infinity
                ia = abs(ia);
            else
                error('Invalid input!')
            end
            latps.ks(ia) = 0; % ignore momentum along ia-th dim
            
            [Ra, ind] = sort(Rs(:,ia));
            step = find(Ra(2:end)-Ra(1:end-1)>0.01);
            rdist = diff([0;step;nR]);
            cRs = mat2cell(Rs(ind,:),rdist); % each cell contains the same Rs(:,ia)
            maxRa = abs(cRs{1}(1,ia));
            
            latps.ns(ia) = 2*maxRa;
            latps.frac = true;
            [ham, pos, part, dof] = h.onLattice(latps);
            
            ind1 = cat(1,part{pos(:,ia)>=maxRa});
            ind0 = find(pos(:,ia)<maxRa);
            pos = pos(ind0,:);
            part = part(ind0);
            ind0 = cat(1,part{:});
            if 2*length(ind0)~=size(ham,1)
                error('Error extracting site indices of principal layer!')
            end       

            h0 = ham(ind0,ind0);
            if norm(ham(ind1,ind0)-ham(ind0,ind1)','fro')/nnz(ham(ind1,ind0))<1e-10 % hermitian case
                if s>0
                    v = ham(ind1,ind0);
                else
                    v = ham(ind0,ind1);
                end
            else % non-hermitian case
                if s>0
                    v = {ham(ind1,ind0), ham(ind0,ind1)};
                else
                    v = {ham(ind0,ind1), ham(ind1,ind0)};
                end
            end
            dof = dof(ind0);
        end

    end
    
    methods (Access=private, Hidden=true)
        function hn = extract(h, inds)
            hn = TBHamiltonian(h.Vectors(inds,:), h.Matrices(:,inds), ...
                               h.DoF, h.UnitCell);
        end
    
        function ham = ham1D(h, n, k, pbc, Bz)
        % Construct Hamiltonian of 1D system along the first dim in h.
        % Note the output is not necessarily hermitian because the 1D
        % system can have a finite displacement in other dimensions.
        % n is the number of unit cells in the 1D systems
        % k (in units of 2pi/a) is momentum (default 0; disabled if Bz nonempty)
        % pbc (default 0) sets periodical boundary condition by being multiplied to hopping terms
        % Bz is in units of Tesla. When Bz is nonempty, use uniform
        % Landau gauge A=(0,Bz*x,0...), and the first vector in the unit
        % cell basis is required to be parallel to (1,0,0...) and all other
        % vectors have to be perpendicular to the first one to allow
        % possible k in the other dimensions.
            if nargin<5
                Bz = [];
            end
            if nargin<4
                pbc = 0;
            else
                pbc = pbc(1);
            end
            if nargin<3
                k = 0;
            else
                k = k(1);
            end
            n = n(1);
            if n<=0
                error('Invalid size of lattice!')
            end
            
            R = h.Vectors;
            H = h.Matrices;
            nR = h.NVecs;
            dH = h.NDMat;
            
            if size(R,2)>1 && nnz(R(:,2:end)-repmat(R(1,2:end),nR,1))>0
                error('Input Rs must differ only in the first dimension!')
            end
            
            ham = sparse(n*dH,n*dH);
            if isempty(Bz)
                for ii = 1:nR
                    Rii = R(ii,1);
                    Hii = exp(-2i*pi*Rii*k)*reshape(H(:,ii),dH,dH);
                    ham = ham + kron(spdiags(ones(n,1),-Rii,n,n),Hii);
                    if pbc*Rii~=0
                        if Rii<0
                            p = pbc;
                        else
                            p = conj(pbc);
                        end
                        if abs(Rii)>=n
                            mRii = mod(Rii,n);
                            if mRii==0
                                ham = ham + p*kron(spdiags(ones(n,1),0,n,n),Hii);
                            else
                                ham = ham + p*kron(spdiags(ones(n,1),-mRii,n,n),Hii);
                                ham = ham + p*kron(spdiags(ones(n,1),n-mRii,n,n),Hii);
                            end
                        else
                            ham = ham + p*kron(spdiags(ones(n,1),sign(Rii)*n-Rii,n,n),Hii);
                        end
                    end
                end
            else % k will not be used because translational symmetry is broken
                eoverh = Utility.PhysConst.e/(Utility.PhysConst.h*1e20);
                uc = h.UnitCell;
                if ~isempty(uc)
                    bas = uc.Basis;
                    if nnz(bas(1,:))>1 || nnz(bas(:,1))>1 || ~bas(1,1)
                        error('First vector in unit cell basis has to be along x and perpendicular to other vectors.')
                    end
                    pos = uc.Sites*bas;
                    nst = uc.NSites;
                    lp = cellfun(@length,h.Partition);
                    mx0 = repelem(reshape(repmat(pos(:,1),nst,1)+repelem(pos(:,1),nst,1), ...
                                          nst,nst)/2, ...
                                  lp,lp); % mean value of x
                    dy0 = repelem(reshape(repmat(pos(:,2),nst,1)-repelem(pos(:,2),nst,1), ...
                                          nst,nst), ...
                                  lp,lp);
                else
                    bas = eye(h.NDims);
                    mx0 = 0;
                    dy0 = 0;
                end
                for ii = 1:nR
                    Rii = R(ii,1);
                    Hii = reshape(H(:,ii),dH,dH);
                    dR = R(ii,:)*bas;
                    dy = dy0+dR(2);
                    for ix = 1:n
                        mx = mx0+dR(1)/2+(ix-(1+n)/2)*bas(1,1);
                        if isscalar(mx)
                            Hiin = exp(2i*pi*Bz*mx*dy*eoverh)*Hii;
                        else
                            Hiin = exp(2i*pi*Bz*mx.*dy*eoverh).*Hii;
                        end
                        if ix+Rii>0 && ix+Rii<=n
                            ham = ham + kron(sparse(ix+Rii,ix,1,n,n),Hiin);
                        elseif pbc~=0
                            if mod(Bz*n*bas(1,1)*dR(2)*eoverh,1)>1e-10
                                warning('Bz and PBC are not compatible!')
                            end
                            if Rii<0
                                ham = ham + pbc*kron(sparse(mod(ix+Rii-1,n)+1,ix,1,n,n),Hiin);
                            else
                                ham = ham + conj(pbc)*kron(sparse(mod(ix+Rii-1,n)+1,ix,1,n,n),Hiin);
                            end
                        end
                    end
                end
            end
        end

    end
    
    methods (Static)
        function [t, basis] = sk(bond, dr)
        % Return the Slater-Koster coefficient matrix in real basis.
        % 
        % Parameters
        % ----------
        % j1, j2 : integers
        %     orbital angular momentum quantum numbers
        % m : integer
        %     component of angular momentum around axis dr (sigma, pi, delta etc.)
        % dr : 3-vector like
        %     displacement vector from 2 to 1
        % 
        % Notes
        % -----
        % The SK coefficient matrix is real and only depends on the angles of dr.
        % It has the following symmetry:
        % sk(j1, j2, m, dr) == sk(j2, j1, m, dr).T
        %                   == (-1)**(j1+j2) sk(j1, j2, m, -dr)
        %                   == (-1)**(j1+j2) sk(j2, j1, m, -dr).T
        % The last equality combined with the Hamiltonian hermicity implies that,
        % the two-center integrals, which only depends on |dr|, satisfy
        % (j2,j1,m) = (-1)**(j1+j2) conj((j1,j2,m))
        % Intuitively this is a consequence of different parities of j1 and j2.
        % 
        % Reference
        % ---------
        % J. C. Slater and G. F. Koster, Phys. Rev. 94, 1498 (1954).
        % R. R. Sharma, Phys. Rev. B 19, 2813 (1979).
            if isscalar(bond) && (nargin<2 || norm(dr)==0) % onsite term
                j1 = Utility.iif(ischar(bond),Utility.orbitmap(bond),bond);
                t = eye(2*j1+1);
                [~, bas] = Symmetry.sph2cub(j1);
                basis = {bas(1,:);bas(1,:)};
                return
            end
            if ischar(bond)
                try
                    j1 = Utility.orbitmap(bond(1));
                    j2 = Utility.orbitmap(bond(2));
                    m = Utility.bondmap(bond(3:end));
                catch
                    error('Invalid input string for bond!')
                end
            else
                j1 = bond(1);
                j2 = bond(2);
                m = bond(3);
            end
            if any([j1, j2, m] < 0) || any(mod([j1, j2, m], 1))
                error('j and m must be non-negative integers')
            end
            if length(dr)~= 3 || norm(dr)==0
                error('dr must be a nonzero three dimensional vector')
            end

            t = zeros(2*j1+1, 2*j2+1);
            if m > min(j1, j2), return, end

            phi = atan2(dr(2), dr(1));
            theta = atan2(norm(dr(1:2)), dr(3));

            for j3 = abs(j1-j2):(j1+j2)
                if mod(j1+j2+j3, 2) % plus/minus m cancel out
                    continue
                elseif m == 0
                    mult = 1;
                else
                    mult = 2; % doubling from plus/minus m
                end
                w1 = Symmetry.wigner3j([j1, j2, j3], [-m, m, 0]);
                if w1 == 0, continue, end
                Ys = Symmetry.sphharm(j3, phi, theta);
                a = zeros(2*j1+1, 2*j2+1);
                for i1 = 1:2*j1+1
                    m1 = j1+1-i1;
                    for i2 = 1:2*j2+1
                        m2 = j2+1-i2;
                        if abs(m2-m1) > j3
                            continue
                        elseif m2-m1 < 0
                            y = (-1)^(m2-m1) * conj(Ys(m1-m2+1));
                        else
                            y = Ys(m2-m1+1);
                        end
                        w2 = Symmetry.wigner3j([j1, j2, j3], [-m1, m2, m1-m2]);
                        a(i1, i2) = (-1)^m1 * w2 * y;
                    end
                end
                t = t + a * sqrt(4*pi*(2*j3+1)) * w1 * mult;
            end

            % transform to real (cubic harmonic) basis
            [U1, bas1] = Symmetry.sph2cub(j1);
            [U2, bas2] = Symmetry.sph2cub(j2);
            t = U1 * t * U2' * (-1)^m;
            if any(reshape(abs(imag(t)),[],1) > 1e5*eps)
                warning(['max imaginary part in SK coefficients is ',...
                      num2str(max(reshape(abs(imag(t)),[],1)))])
            end
            t = real(t);
            basis = {bas1(1,:);bas2(1,:)};
        end
    
        function [Rs, HR] = hamSK(uc, dof, terms)
        % Construct Slater-Koster TB Hamiltonian.
        % terms is cell array with each row representing one term in format:
        % {site label or {pair of labels in cell},
        %  neighbor order,
        %  string (e.g. 'pppi') or [j1,j2,m] or function handle,
        %  coupling strength}
        %   if the first element is one single site label, then only construct the
        %   terms from that site to its corresponding neighbors (not necessarily hermitian),
        %   in this case j2 corresponds to that site;
        %   if the first element is a pair of labels, then construct the
        %   terms connecting the two sites (hermitian) after checking they 
        %   are consistent with the order of neighbors, in this case j1(j2) 
        %   corresponds to the first(second) site. Note if both j1,j2 and
        %   orbitals are the same, this form can lead to double counting.
            Rs = [];
            HR = [];
            nd = length(dof);
            for ii = 1:size(terms,1)
                t = terms(ii,:);
                if iscell(t{1})
                    paired_sites = true;
                    [lab1, lab2] = t{1}{:};
                else
                    paired_sites = false;
                    lab2 = t{1};
                end
                neighbor = t{2};
                [ucrs, labs, drs] = uc.findNeighbors(lab2,neighbor);
                for inb = 1:size(ucrs,1)
                    if paired_sites && ~isequal(labs{inb},lab1)
                        continue % ignore term with inconsistent j1 and j2
                    end
                    bond = t{3};
                    if isa(bond,'function_handle')
                        [tsk, bas] = bond(drs(inb,:));
                    else
                        [tsk, bas] = TBHamiltonian.sk(bond, drs(inb,:));
                    end
                    ind1 = TBHamiltonian.findIndices(dof,...
                        'Site',labs{inb},'Orbital',bas{1});
                    ind2 = TBHamiltonian.findIndices(dof,...
                        'Site',lab2,'Orbital',bas{2});
                    HRi = zeros(nd);
                    HRi(ind1,ind2) = t{4}*kron(tsk,eye(length(ind1)/size(tsk,1)));
                    HR = cat(3,HR,HRi);
                    Rs = [Rs; ucrs(inb,:)]; %#ok<*AGROW>
                    if paired_sites && neighbor>0
                        HR = cat(3,HR,HRi');
                        Rs = [Rs; -ucrs(inb,:)];
                    end
                end
            end
        end
        
        function [Rs, HR] = importDFT(filename, nhead, hasDegR)
        % Import TB model from DFT
            if nargin<3
                hasDegR = true;
            end
            
            delimiter = ' ';
            formatSpec = '%f';
            fileID = fopen(filename,'r');

            NB = cell2mat(textscan(fileID, formatSpec, 1, ...
                        'HeaderLines', nhead, ...
                        'Delimiter', delimiter, ...
                        'MultipleDelimsAsOne', true, ...
                        'ReturnOnError', false));
            NR = cell2mat(textscan(fileID, formatSpec, 1, ...
                        'Delimiter', delimiter, ...
                        'MultipleDelimsAsOne', true, ...
                        'ReturnOnError', false));
            if hasDegR
            DegR = cell2mat(textscan(fileID, formatSpec, NR, ...
                        'Delimiter', delimiter, ...
                        'MultipleDelimsAsOne', true, ...
                        'ReturnOnError', false));
            end
            data = cell2mat(textscan(fileID, formatSpec, ...
                        'Delimiter', delimiter, ...
                        'MultipleDelimsAsOne', true, ...
                        'ReturnOnError', false));

            fclose(fileID);

            data = reshape(data,7,[]);
            if size(data,2)~=NR*NB^2
                Rs0 = data(1:3,:).';
                iRs0 = [find(sum(abs(diff(Rs0)),2)); size(Rs0,1)];
                Rs = Rs0(iRs0,:);
                NR = size(Rs,1);
                ham = @(m) reshape(sparse(m(1,:),m(2,:),m(3,:)+1i*m(4,:),NB,NB),[],1);
                HR = cellfun(ham, mat2cell(data(4:end,:),4,[iRs0(1);diff(iRs0)]),...
                             'UniformOutput', false);
                HR = cell2mat(HR);
            else
                Rs = reshape(permute(reshape(data(1:3,:),3,NB^2,NR),[1,3,2]),3*NR,NB^2);
                if norm(Rs-kron(Rs(:,1),ones(1,NB^2)),'fro')>1e-8
                    error('Unexpected format!')
                end
                Rs = reshape(Rs(:,1),3,NR).';
                HR = reshape(data(end-1,:)+1i*data(end,:),NB^2,NR);
            end

            if hasDegR
                HR = HR./kron(reshape(DegR(1:NR),1,[]),ones(NB^2,1));
            end
        end

        function [Rs, HR] = combineSameR(Rs0, HR0, cutoff)
        % Sort (from last to first dim) Rs0 and combine HR0 for same R
            if nargin<3
                cutoff = 1e-10;
            end
            nd = size(Rs0,2);
            nR = size(Rs0,1);
            [Rs0, indR] = sortrows(Rs0,nd:-1:1);
            HR0 = sparse(reshape(HR0,[],nR));
            HR0 = HR0(:,indR);
            step = [find(sum(abs(diff(Rs0,1,1)),2)>0.01);nR]; % last indices before R change
            if length(step)==nR
                Rs = Rs0;
                HR = HR0;
            else
                lens = diff([0;step]); % length of same R
                for ii = find(lens>1).'
                    HR0(:,step(ii)) = sum(HR0(:,step(ii)+(1-lens(ii):0)),2); %#ok<SPRIX>
                end
                Rs = Rs0(step,:);
                HR = HR0(:,step);
            end
            nzHR = find(sum(abs(HR),1)>cutoff);
            if length(nzHR)~=length(step)
                Rs = Rs(nzHR,:);
                HR = HR(:,nzHR);
            end
        end
        
        function tf = isHermitian(Rs, HR, tol)
            if nargin<3
                tol = 1e-10;
            end
            [Rs, HR] = TBHamiltonian.combineSameR(Rs, HR);
            nR = size(Rs,1);
            dH = sqrt(size(HR,1));
            for iR = 1:nR
                R = Rs(iR,:);
                iR2 = Utility.findvec(-R,Rs);
                if isempty(iR2)
                    tf = false;
                    return
                end
                H1 = reshape(HR(:,iR),dH,dH);
                H2 = reshape(HR(:,iR2),dH,dH);
                if norm(H1-H2','fro')>tol
                    tf = false;
                    return
                end
            end
            tf = true;
        end
        
        function poss = allPosition(pos, part)
        % Obtain positions for all indices from positions of all sites
            n = size(pos,1);
            assert(n==length(part), 'Invalid input size!')
            poss = cell2mat(cellfun(@(p,c)repmat(p,length(c),1),...
                                    num2cell(pos,2), part,...
                                    'UniformOutput', false));
            poss(cell2mat(part),:) = poss;
        end
        
        function [dof, part] = convertDoF(dofin)
        % Convert dof from input format to stored format.
        % Output dof is a struct array: dof.(name)=(value)
        % Examples of input format:
        % dof = {...
        %     {'P1',1:4};...
        %     {'P2',5:8}};
        % dof = {...
        %     {'',1:8,...
        %         'Spin',{'u','d'},...
        %         'Orbital',{'s','px','py','pz'}}};
        % dof = {...
        %     {'P1',[1:4,17:20],...
        %         'Nambu',{'e','h'},...
        %         'Spin',{'u','d'},...
        %         'Orbital',{'s','pz'}};...
        %     {'P2',[5:8,21:24],...
        %         'Nambu',{'e','h'},...
        %         'Spin',{'u','d'},...
        %         'Orbital',{'s','pz'}};...
        %     {'P3',[9:12,25:28],...
        %         'Nambu',{'e','h'},...
        %         'Spin',{'u','d'},...
        %         'Orbital',{'s','pz'}};...
        %     {'P4',[13:16,29:32],...
        %         'Nambu',{'e','h'},...
        %         'Spin',{'u','d'},...
        %         'Orbital',{'s','pz'}}};
        % dof = {...
        %     {'P1',[1:4,17:20],...
        %         'Nambu',{'e','h'},...
        %         'Spin',{'u','d'},...
        %         'Orbital',{'s','pz'}};...
        %     {'P2',[5:8,21:24]};...
        %     {'P3',[9:12,25:28]};...
        %     {'P4',[13:16,29:32]}};
        % dof = {...
        %     {'O',[1:8,17:24],...
        %         'Nambu',{'e','h'},...
        %         'Spin',{'u','d'},...
        %         'Orbital',{'s','px','py','pz'}};...
        %     {'P',[9:14,25:30],...
        %         'Nambu',{'e','h'},...
        %         'Spin',{'u','d'},...
        %         'Orbital',{'s','px','py'}};...
        %     {'Q',[15:16,31:32],...
        %         'Nambu',{'e','h'},...
        %         'Spin',{'u','d'},...
        %         'Orbital',{'s'}}};
            for ii = 1:length(dofin)
                dofii = dofin{ii};
                if length(dofii)>2 || ii==1
                    headings = {'Site',dofii{3:2:end}};
                    labels = {dofii(1),dofii{4:2:end}};
                    ndofs = cellfun(@(c)length(c),labels);
                    labels = arrayfun(@(c,i)repmat(repelem(reshape(c{:},[],1),...
                                                         prod(ndofs(i+1:end)),1),...
                                                 prod(ndofs(1:i-1)),1),...
                                    labels,1:length(labels),'UniformOutput',false);
                    labels = cat(2,labels{:});
                else
                    labels(:,1) = dofii(1);
                end
                dof(dofii{2}) = cell2struct(labels,headings,2); 
            end
            dof = reshape(dof,[],1);
            part = cellfun(@(c)reshape(c{2},[],1),dofin,'UniformOutput',false);
        end
        
        function inds = findIndices(dof, varargin)
        % Find indices according to the field and value pairs of dof.
        % Input values can be one value (string or numeric matrix) or
        % cell array of values
        % Eg. findIndices(dof, ...
        %                 'field1_name', 'filed1_value',...
        %                 'field2_name', {filed2_value1,filed2_value2})
            if mod(nargin-1,2)~=0
                error('Invalid number of input!')
            end
            if ~isstruct(dof)
                try
                    [dof, ~] = TBHamiltonian.convertDoF(dof);
                catch
                    error('Invalid input dof!')
                end
            end
            inds = (1:length(dof)).';
            nf = (nargin-1)/2;
            ii = 1;
            while ii<=nf && ~isempty(inds)
                l = varargin{2*ii-1};
                v = varargin{2*ii};
                ind = @(v1) reshape(find(cellfun(@(v2)isequal(v1,v2),{dof.(l)})),[],1);
                if ~iscell(v)
                    indii = ind(v);
                else
                    indii = cell2mat(cellfun(@(v1)ind(v1),reshape(v,[],1),'UniformOutput',false));
                end
                inds = inds(indii);
                dof = dof(indii);
                ii = ii+1;
            end
        end
        
        function op = fullOperator(dof, op0)
        % Eg. fullOperator(dof,...
        %                  {'Spin',  {'u','d'}, [0,1;1,0];...
        %                   'Nambu', {'e','h'}, [1,0;0,0]})
            assert(isstruct(dof), 'dof must be a struct!')
            ndof = length(dof);
            if nargin==1
                op = speye(ndof);
                return
            end
            fns = op0(:,1);
            fns2 = setdiff(fieldnames(dof),fns);
            dofop = cell2struct(Utility.cartprod(op0{:,2}),fns,2);
            indall = 1:ndof;
            op = zeros(ndof);
            while ~isempty(indall)
                dofall = dof(indall);
                d0 = rmfield(dofall(1),fns);
                ind0 = Utility.findstruct(d0, rmfield(dofall,fns));
                if length(ind0)~=length(dofop)
                    error('Incompatible dof!')
                end
                ord = Utility.findstruct(dofop, rmfield(dofall(ind0),fns2));
                ind = indall(ind0(ord));
                op(ind,ind) = Utility.kronall(op0{:,3});
                indall = setdiff(indall,ind);
            end
            op = sparse(op);
        end
        
        function ham = addOnsiteTerm(ost, ham, pos, part, dof)
        % ost(pos,dof) returns the on-site term matrix
        % or, ost(pos) returns a function handle taking (pos,dof) as input
        % the later form is to allow different functions for different pos
            n = size(pos,1);
            if nargin(ost)==1 % passing function handle to delay conditional evaluation
                ostf = @(f,p,d) f(p,d);
                osts = cellfun(ostf,...
                        arrayfun(@(ii)ost(pos(ii,:)),1:n,'UniformOutput',false),...
                        num2cell(pos,2),...
                        arrayfun(@(ii)dof(part{ii}),1:n,'UniformOutput',false),...
                        'UniformOutput',false);
            elseif nargin(ost)==2
                osts = arrayfun(@(ii)ost(pos(ii,:),dof(part{ii})),...
                        1:n,'UniformOutput',false);
            else
                error('Unrecoganized input ost!')
            end
            ind = cat(1,part{:});
            ham(ind,ind) = ham(ind,ind)+blkdiag(osts{:});
        end
        
        function ham = rescaleSiteCoupling(r, ham, part, mask1, mask2)
        % to rescale coupling matrix by a factor/function r
            part1 = part(mask1);
            ind1 = cat(1,part1{:});
            part2 = part(mask2);
            ind2 = cat(1,part2{:});
            if isa(r,'function_handle')
                ham(ind1,ind2) = r(ham(ind1,ind2));
            else
                ham(ind1,ind2) = r*ham(ind1,ind2);
            end
            ham(ind2,ind1) = ham(ind1,ind2)';
        end
        
        function HR = TR1(HR0)
        % time reversal transformation for one spin species
        %   simply complex conjugate is applied
            HR = conj(HR0);
        end
        
        function HR = minusTR2(HR0, dof_or_UT)
        % time reversal transformation for two spin species, multiplied by
        % a minus sign at last, as used in BdG Hamiltonian
        % T = i*sigma_y*K is used for the TR operator if dof is provided,
        % T = UT*K is used if UT is provided
            if isempty(dof_or_UT)
                UT = kron([0,1;-1,0],eye(size(HR0,1)/2));
            elseif isstruct(dof_or_UT)
                dof = dof_or_UT;
                UT = TBHamiltonian.fullOperator(dof, {'Spin',{'u','d'},[0,1;-1,0]});
            elseif isnumeric(dof_or_UT) && isequal(size(dof_or_UT),size(HR0))
                UT = dof_or_UT;
            else
                error('Invalid input dof_or_UT!')
            end
            HR = UT*conj(HR0)*UT;
        end
    end
end