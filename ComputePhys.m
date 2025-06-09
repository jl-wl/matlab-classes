classdef ComputePhys
% ComputePhys contains routines to compute physical quantities.
    
    properties
    end
    
    methods (Static)
        function [qs, corr, rmask] = rhoCorr(ks, rhos, mask, extbz)
        % Correlation function of spectral (and spin) density in k-space.
        % Inputs:
        %   ks (cell or numeric array): momenta for rhos;
        %   rhos (numeric array): spectral density, or spin and spectral
        %     densities (arranged as [s1,s2,s3,s0] on the last dimension)
        %     at a constant energy;
        %   mask (struct): mask setting for rhos;
        %   extbz (numeric array): extended BZ (boundary values),
        %     here BZ is extended by padding with zeros.
        % Outputs:
        %   qs (cell or numeric array): displacement momenta for corr;
        %   corr (numeric array): correlation function,
        %     if rhos contains only spectral density, then corr is jdos;
        %     if rhos contains both spin and spectral densities, then corr
        %     is cat(nd+1,jdos,ssp), ssp being spin scattering probability;
        %   rmask (logic or numeric array): mask weights for rhos
            if iscell(ks)
                nd = length(ks); % spacial dimension
            elseif isnumeric(ks)
                nd = 1;
                ks = {ks};
            else
                error('Invalid input ks!')
            end
            kpts = reshape(cellfun(@length,ks),1,[]);
            
            nrho = size(rhos);
            ndrho = length(nrho);
            if (ndrho~=nd && ndrho~=nd+1) || ~isequal(kpts,nrho(1:nd))
                error('Size of ks does not match size of rhos!');
            end
            rhos = num2cell(rhos,1:nd);
            rho0 = rhos{end}; % use the last slice for spectral density
            
            if nd==1
                kpts = [kpts,1];
            end

            if nargin<4
                qpts = kpts;
                qs = ks;
            else
                qpts = cellfun(@(b,k)ceil(2*b/(k(2)-k(1)))+1,num2cell(extbz),ks);
                if nd==1
                    qs = linspace(-extbz,extbz,qpts);
                    qpts = [qpts,1];
                else
                    qs = arrayfun(@(b,n)linspace(-b,b,n),extbz,qpts,...
                                  'UniformOutput',false);
                end
            end
            
            if nargin<3 || isempty(mask)
                mask.type = 'none';
            end
            switch mask.type
                case 'none'
                    rmask = ones(kpts);
                case 'min' % use a lower bound for spectral density
                    rmask = rho0>=mask.value;
                case 'max' % use an upper bound for spectral density
                    rmask = rho0<=mask.value;
                case 'region' % select union of specific regions in k-space
                    % each row of mask.value sets a region by [k1min, k1max, k2min, k2max, ...]
                    nreg = size(mask.value,1); % number of regions
                    rmask = zeros(kpts);
                    for ireg = 1:nreg
                        ms = cellfun(@(k,bs)(k>=bs(1))&(k<=bs(2)),...
                                     ks,mat2cell(mask.value(ireg,:),1,2*ones(1,nd)),...
                                     'UniformOutput',false); % masks for each dimension
                        cms = cell(nd,1);
                        [cms{:}] = ndgrid(ms{:});
                        rmask = rmask | all(cat(nd+1,cms{:}),nd+1);
                    end
                case 'regionminmax' % combination of 'region', 'min', and 'max'
                    if isfield(mask,'regvalue')
                        nreg = size(mask.regvalue,1); % number of regions
                        rmask = zeros(kpts);
                        for ireg = 1:nreg
                            ms = cellfun(@(k,bs)(k>=bs(1))&(k<=bs(2)),...
                                         ks,mat2cell(mask.regvalue(ireg,:),1,2*ones(1,nd)),...
                                         'UniformOutput',false); % masks for each dimension
                            cms = cell(nd,1);
                            [cms{:}] = ndgrid(ms{:});
                            rmask = rmask | all(cat(nd+1,cms{:}),nd+1);
                        end
                    else
                        rmask = ones(kpts);
                    end
                    if isfield(mask,'minvalue')
                        rmask = rmask & (rho0>=mask.minvalue);
                    end
                    if isfield(mask,'maxvalue')
                        rmask = rmask & (rho0<=mask.maxvalue);
                    end
                case 'rescale' % rescale weight in region by factor of w
                    % weight outside the specified regions is 1;
                    % w is stored in last entry of each row of mask.value
                    nreg = size(mask.value,1); % number of regions
                    rmask = ones(kpts);
                    for ireg = 1:nreg
                        w = mask.value(ireg,end); % scaling factor
                        ms = cellfun(@(k,bs)(k>=bs(1))&(k<=bs(2)),...
                                     ks,mat2cell(mask.value(ireg,1:end-1),1,2*ones(1,nd)),...
                                     'UniformOutput',false); % masks for each dimension
                        cms = cell(nd,1);
                        [cms{:}] = ndgrid(ms{:});
                        rmaski = ones(kpts)-(1-w)*all(cat(nd+1,cms{:}),nd+1);
                        rmask = rmask.*rmaski;
                    end
                case 'envelop' % rescale by an envelop function
                    fenv = mask.func; % fenv(k1,k2,...) gives scaling factor
                    gks = cell(nd,1);
                    [gks{:}] = ndgrid(ks{:});
                    rmask = fenv(gks{:});
            end

            jall = cellfun(@jrho,rhos,'UniformOutput',false);
            if isscalar(jall)
                corr = jall{1};
            else
                ssp = 0.5*sum(cat(nd+1,jall{:}),nd+1);
                corr = cat(nd+1,jall{end},ssp);
            end

            function j = jrho(rho)
                rho = rho.*rmask;
                j = Utility.chop(fftshift(fftn(fftn(rho,qpts).*ifftn(rho,qpts)))/prod(qpts));
                % note that this method assumes rhos is periodical in qs
            end
        end
    
        function a = spf(eni, ens, wfs, se)
        % Spectral function (diagonal) at real(eni) from eigenstates
        % (ens, wfs) with optional self-energy
            assert(imag(eni)>0,'Imaginary part of eni must > 0!')
            if max(real(eni))+imag(eni)>max(ens) || min(real(eni))-imag(eni)<min(ens)
                warning('Input energy range of eigenstates not large enough!')
            end
            inve = Utility.spd(1./(eni-reshape(ens,[],1)),0);
            if nargin>3
                [r,c] = find(se);
                ind1 = union(r,c);
                ind2 = setdiff((1:size(wfs,1)).',ind1);
                se = se(ind1,ind1);
                wfs_inve = wfs*inve;
                gf11 = wfs_inve(ind1,:)*wfs(ind1,:)';
                gf21 = wfs_inve(ind2,:)*wfs(ind1,:)';
                gf12 = wfs_inve(ind1,:)*wfs(ind2,:)';
                gf22d = sum(wfs_inve(ind2,:).*conj(wfs(ind2,:)),2);
                gf11 = gf11/(eye(length(ind1))-se*gf11);
                gf21 = gf21*(eye(length(ind1))+se*gf11);
                gf22d = gf22d+sum((gf21*se).*gf12.',2);
                a(ind1,1) = -imag(diag(gf11))/pi;
                a(ind2,1) = -imag(gf22d)/pi;
            else
                a = -imag(sum((wfs*inve).*conj(wfs),2))/pi;
            end
            if nnz(a<0)
                warning('Negative spectral weight found!')
            end
        end
        
        function rho = ldos(es, pos, es0, pos0, rho0, fwe, fwp)
        % Local density of states in measurement.
        % pos0 and rho0 can be cell array in case of multiple parts;
        % size of raw rho0 data is size(pos0,1)*length(es0);
        % fwe is function (or cell array of functions) for weight as a 
        % function of energy difference; fwp is function (or cell array of
        % functions) for weight as a function of measured position (pos)
        % and initial position (pos0)
            ne = length(es);
            np = size(pos,1);
            de = repmat(reshape(es,1,[]),length(es0),1) - repmat(reshape(es0,[],1),1,ne);
            if iscell(pos0)
                n = length(rho0); % number of parts
            else
                n = 1;
                pos0 = {pos0};
                rho0 = {rho0};
            end
            rho = zeros(np,ne);
            for ii = 1:n
                p0 = pos0{ii};
                np0 = size(p0,1);
                p = repmat(pos,np0,1);
                p0 = repelem(p0,np,1);
                if iscell(fwe)
                    we = fwe{ii};
                else
                    we = fwe;
                end
                if iscell(fwp)
                    wp = fwp{ii};
                else
                    wp = fwp;
                end
                rho = rho + reshape(wp(p,p0),np,np0)*rho0{ii}*we(de);
            end
        end
        
        function varargout = ev(ops, wfs)
        % Expectation values of operators
            s = size(wfs);
            if isnumeric(ops)
                ops = {ops};
            end
            if ~iscell(ops) || ~all(cellfun(@isnumeric,ops))
                error('Invalid input for operators!')
            end
            nop = length(ops);
            vals = cell(1,nop);
            for ii = 1:nop
                op = ops{ii};
                dop = size(op,1);
                if size(op,2)~=dop || mod(s(1),dop)~=0
                    error('Invalid operator!')
                end
                wfs = reshape(wfs,dop,[]);
                vals{ii} = reshape(Utility.chop(sum(conj(wfs).*(op*wfs),1)),[s(1)/dop,s(2:end)]);
            end
            if nargout==1 && nop>1
                varargout{1} = vals;
            else
                for ii = 1:nargout
                    varargout{ii} = vals{ii}; 
                end
            end
        end

        function [phis, ens] = wn(fham, VC, ks)
        % Winding number for Hamiltonian with chiral symmetry.
        %   fham(k) constructs the Hamiltonian;
        %   VC is/constructs transformation that diagonalize chiral op;
        %   ks is list of k.
            [phis, ens] = arrayfun(@p,ks,'UniformOutput',false);
            phis = cell2mat(reshape(phis,1,[]));
            ens = cell2mat(reshape(ens,1,[]));
            
            function [phi, es] = p(k)
                if isa(VC,'function_handle')
                    V = VC(k);
                else
                    V = VC;
                end
                hamk = V'*fham(k)*V;
                d = size(hamk,1)/2;
                if any(Utility.chop(hamk(1:d,1:d)),'all') ...
                        || any(Utility.chop(hamk(d+1:end,d+1:end)),'all')
                    error('Diagonal block is nonzero!')
                end
                [V1,S,V2] = svd(hamk(1:d,d+1:2*d));
                phi = sort(angle(eig(V1*V2')));
                es = sort(diag(S));
            end
        end

        function [phis, gaps] = wl(fham, nocc, ks)
        % Wilson loop spectrum and gaps between occupied and unoccupied bands.
        %   fham(k) constructs the Hamiltonian;
        %   nocc is the number of occupied states;
        %   ks can be a COLUMN list of k, or cell array of such lists.
            if ~iscell(ks)
                [phis, gaps] = w(ks);
            else
                [phis, gaps] = cellfun(@(kl)w(kl),ks,'UniformOutput',false);
                phis = reshape(cell2mat(phis),nocc,[]);
            end
            
            function [phi, gap] = w(kl)
%                 tic
                nkl = size(kl,1);
                gap = zeros(nkl,1);
                for ikl = 1:nkl
                    [wfs, ens] = Utility.eigen(full(fham(kl(ikl,:))));
                    gap(ikl) = ens(nocc+1)-ens(nocc);
                    wfs = wfs(:,1:nocc);
                    if ikl==1
                        wfs1 = wfs;
                        p = wfs';
                    else
                        p = p*wfs*wfs'; %#ok<MHERM>
                    end
                end
                p = p*wfs1;
                phi = angle(eig(p))/pi;
%                 toc
            end
        end
        
        function [bi, ens] = bi2d(ham, poss, nocc)
        % Bott index in 2D %TODO: test
        %   poss is fractional coordinates for all indices
        %   nocc is number of occupied states
            [evecs, ens] = Utility.eigen(full(ham));
            evecs = evecs(:,1:nocc);
            nx = floor(max(poss(:,1)))-floor(min(poss(:,1)))+1;
            ny = floor(max(poss(:,2)))-floor(min(poss(:,2)))+1;
            Vx = evecs'*Utility.spd(exp(2i*pi*poss(:,1)/nx))*evecs;
            Vy = evecs'*Utility.spd(exp(2i*pi*poss(:,2)/ny))*evecs;
            bi = imag(sum(log(eig(Vy*Vx*Vy'*Vx'))))/(2*pi);
        end
        
        function tm = tMat(emh, v)
        % emh = en+i*eta-h0
        % v is from i-th layer to (i+1)-th
            d = size(emh,1);
            t = emh\v; tm = t;
            tp = emh\v'; tpp = tp;
            while abs(norm(full(t)))+abs(norm(full(tp)))>1e-10
                a = speye(d)-t*tp-tp*t;
                t = a\t*t;
                tp = a\tp*tp;
                tm = tm+tpp*t;
                tpp = tpp*tp;
            end
        end

        function [gk, frk, ks] = sGFlead(nsec, en, t, mu)
        % To compute surface (retarded) Green function for (spinless)
        % normal metal lead.
        % nsec is real-space size (in row vector) for the section of lead.
        % gk is sGF for each k-point corresponding to nsec.
        % frk is function handle for the transformation from k to r,
        % such that g(r1,r2) = fr1k.gk.fr2k'.
        % frk can take multiple r and multiple k (each row is one point).
            if nargin<4
                mu = 0;
            end
            if nargin<3
                t = 1;
            end
            if isreal(en)
                en = en+1i*1e-5;
            end
            en = (en+mu)/(2*t)-1; % energy in normalized form (i.e. -cosine)

            ks = Utility.ndgridVecs(ones(size(nsec)),nsec);
            ens = en - sum(1-cos(ks*diag(pi./(nsec+1))),2);
            gk = exp(-1i*acos(ens))/t; % retarded GF has negative imaginary part

            % orthonormal wf in 1D:
            % \psi_{i_r, i_k} = \sqrt{\frac{2}{L+1}}\sin(\frac{\pi i_r i_k}{L+1})
            % where i_r and i_k are both from 1 to L
            frk = @(r,k) reshape(prod(sin(repmat(r,size(k,1),1).*repelem(k,size(r,1),1)*diag(pi./(nsec+1))),2),size(r,1),size(k,1))*sqrt(prod(2./(nsec+1)));
        end
        
        function rgf = rGFchain(ps,en,se) % TODO
        % To compute retarded GF at either end AND between two ends;
        % only sites connected to leads are maintained.
        % Use this function only when leads are at 1D ends.

            nx = ps.xLength;
            ny = ps.yLength;
            dof = ps.numComponents;
            db = ny*dof;
            fh = ps.xHamIntraLayer; % Ham at layer l
            fv = ps.xHamInterLayer; % Ham from layer l+1 to l

            if imag(en)<1e-10
                en = real(en) + 1e-10i;
            end

            sel = zeros(db);
            ser = zeros(db);
            inds = zeros(1,0);
            leads = ps.Leads;
            nl = length(leads);
            for il = 1:nl
                l = leads(il);
                ind = reshape(kron(1-dof:0,ones(length(l.index),1))' ...
                    +kron(l.index*dof,ones(dof,1)),1,[]);
                switch l.location(2)
                    case '-'
                        sel(ind,ind) = se{il};
                        inds = [inds ind]; %#ok<*AGROW>
                    case '+'
                        ind = mod(ind-1,db)+1;
                        ser(ind,ind) = se{il};
                        inds = [inds ind+db];
                end
            end

            emh = en*eye(db)-full(fh(ps,1));

            if nx<2
                rgf = inv(emh-sel-ser);
                rgf = rgf(inds,inds);
                return
            end

            rgf = zeros(2*db);
            lt = 1:db;
            rt = db+1:2*db;

            rgf(lt,lt) = inv(emh-sel);
            rgf(rt,rt) = rgf(lt,lt);
            rgf(lt,rt) = rgf(lt,lt);
            rgf(rt,lt) = rgf(lt,lt);

            for ix = 2:nx
                emh = en*eye(db)-full(fh(ps,ix));
                vx = full(fv(ps,ix-1));
                if ix<nx
                    rgf(rt,rt) = inv(emh-vx'*rgf(rt,rt)*vx);
                else
                    rgf(rt,rt) = inv(emh-vx'*rgf(rt,rt)*vx-ser);
                end
                rgf(rt,lt) = rgf(rt,rt)*vx'*rgf(rt,lt);
                rgf(lt,lt) = rgf(lt,lt) + rgf(lt,rt)*vx*rgf(rt,lt);
                rgf(lt,rt) = rgf(lt,rt)*vx*rgf(rt,rt);
            end

            rgf = rgf(inds,inds);
        end
    end
    
end

