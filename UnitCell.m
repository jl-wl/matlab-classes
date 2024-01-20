classdef UnitCell
    %UNITCELL Summary of this class goes here
    %Convention: position/k vector in same row
    
    properties
        Basis % in units of Angstrom
        Sites % fractional coordinates
        Origin % origin in absolute coordinates
        Labels % labels for each site
        NDims
        NSites
    end
    
    methods
        function uc = UnitCell(basis,sites,origin,labels)
        % Note: sites must be fractional coordinates
            if nargin>0
                uc.Basis = basis;
                if abs(det(basis))<1e-12
                    error('Basis is singular!')
                end
            else
                error('Must at least provide basis!')
            end
            nd = size(basis,2);
            if nargin<2
                sites = zeros(1,nd);
            elseif size(sites,2)~=nd
                error('Dimensions mismatch!')
            end
            n = size(sites,1);
            if nargin<3
                origin = zeros(1,nd);
            elseif size(origin,1)~=1 || size(origin,2)~=nd
                error('Dimensions mismatch!')
            end
            if nargin<4 || isempty(labels)
                labels = {};
            elseif length(labels)~=n
                error('Dimensions mismatch!')
            else
                labels = reshape(labels,[n,1]);
            end
            uc.Sites = sites;
            uc.Origin = origin;
            uc.Labels = labels;
            uc.NDims = nd;
            uc.NSites = n;
        end

        function kbas = KBasis(uc)
            kbas = inv(uc.Basis).';
        end
        
        function [ucrs, labs, drs, ds] = findNeighbors(uc, lab0, n, tol)
        % Find n-th (n can be a vector) nearest neighbors.
        % ucrs contains displacement vector in unit cells for each neighbor
        % labs contains label of each neighbor
        % drs / ds contains the absolute displacement vector / distance for each neighbor
        % if n is vector, then ucrs, labs, drs are cell arrays with each
        % cell corresponding to each neighbor order
            if nargin<4
                tol = 1e-8; % tolerance in determining equal distance
            end
            bas = uc.Basis;
            pos = uc.Sites;
            lab = uc.Labels;
            nd = uc.NDims;
            n0 = uc.NSites;
            p0 = pos(strcmp(lab0,lab),:);
            if isempty(p0)
                error([lab0,' cannot be found!'])
            elseif size(lab0,1)>1
                error([lab0,' cannot be uniquely determined!'])
            end
            dp0 = pos - repmat(p0,n0,1);
            nuc0 = ceil(max(n)/n0);
            ds = []; ls = []; ds_old = []; ls_old = [];
            while isempty(ds) || isempty(ds_old) || norm(ds-ds_old)>tol || norm(ls-ls_old)>tol
                ucrs0 = Utility.ndgridVecs(-nuc0*ones(1,nd),nuc0*ones(1,nd));
                nuc = size(ucrs0,1);
                dp = (repmat(dp0,nuc,1)+repelem(ucrs0,n0,1))*bas;
                [dss, ind] = sort(sqrt(sum(dp.^2,2)));
                step = [find(diff(dss)>tol);nuc*n0];
                lens = diff([0;step]);
                indc = mat2cell(ind,lens);
                ucrs = cellfun(@(ii)ucrs0(ceil(ii/n0),:),indc(n+1),'UniformOutput',false);
                labs = cellfun(@(ii)lab(mod(ii-1,n0)+1,:),indc(n+1),'UniformOutput',false);
                drs = cellfun(@(ii)dp(ii,:),indc(n+1),'UniformOutput',false);
                ds_old = ds;
                ls_old = ls;
                ds = dss(step(n+1));
                ls = lens(n+1);
                nuc0 = nuc0+1;
            end
            if length(n)==1
                ucrs = ucrs{1};
                labs = labs{1};
                drs = drs{1};
            end
        end
        
        function ucn = rotateCoordinates(uc0,rot)
            if abs(det(rot))<1e-12
                error('Rotation matrix is singular!')
            end
            ucn = uc0;
            ucn.Basis = uc0.Basis*rot;
            ucn.Origin = uc0.Origin*rot;
        end
        
        function ucn = reconstruct(uc0,new_basis,new_origin)
            bas0 = uc0.Basis;
            pos0 = uc0.Sites;
            org0 = uc0.Origin;
            lab0 = uc0.Labels;
            nd = uc0.NDims;
            n0 = uc0.NSites;
            
            if nargin>1
                basn = new_basis;
            else
                error('Must at least provide new basis!')
            end
            if nargin>2
                orgn = new_origin;
            else
                orgn = zeros(1,nd);
            end
            
            if abs(det(basn))<1e-12
                error('New basis is singular!!')
            end
            if any(nd~=size(basn))
                error('Dimensions mismatch!')
            end
            
            rbasn = basn/bas0; % new basis vectors in old basis
            if norm(rbasn-round(rbasn),'fro')>1e-12
                error('Incompatible basis!')
            end
            rbasn = round(rbasn);

            m = det(rbasn); % ratio of volumn new to old
            if m<0
                warning('New basis has an opposite signature!')
                m = -m;
            end

            % search for translational vectors in old basis that fall into new uc
            mins = sum(rbasn.*(rbasn<0),1); % min in each dim
            maxs = sum(rbasn.*(rbasn>0),1); % max in each dim
            v = Utility.ndgridVecs(mins,maxs,1); % vectors to be searched
            v = v*bas0/basn+1e-8; % convert v to new basis, shift a tiny amount to remove ambiguity
            v = v(all(abs(v-0.5)<0.5,2),:); % keep vectors inside new uc
            if size(v,1)~=m
                error('Error finding all inequivalent points in new uc!')
            end
            
            pos0n = mod((pos0*bas0+repmat(org0-orgn,[n0,1]))/basn,1);
            posn = mod(repmat(pos0n,m,1)+repelem(v,n0,1),1);
            
            if ~isempty(lab0) && m>1
                labn = cat(2,repmat(lab0,m,1),num2cell(repelem((1:m).',n0,1)));
                labn = cellstr(join(string(labn),{'_'}));
            else
                labn = lab0;
            end
            
            ucn = UnitCell(basn,posn,orgn,labn);
        end
        
        function plotSites(uc,ns)
        % ns is vector of number of unit cells in each dimension
            if isempty(uc.Sites) || isempty(uc.Basis)
                error('Sites or basis is not provided!')
            end
            nd = uc.NDims;
            if nargin>1
                if length(ns)~=nd || any(ns<1)
                    error('Invalid input ns!')
                else
                    ns = reshape(ns,1,[]);
                end
            else
                ns = ones(1,nd);
            end
            pos0 = uc.Sites*uc.Basis;
            ucrs = Utility.ndgridVecs(zeros(1,nd),ns-ones(1,nd))*uc.Basis;
            nuc = size(ucrs,1);
            vs = Utility.ndgridVecs(zeros(1,nd),ones(1,nd))*uc.Basis;
            c0 = get(groot,'defaultAxesColorOrder');
            set(groot,'defaultAxesColorOrder',jet(uc.NSites));
            hold on
            switch nd
                case 1
                    for ii = 1:uc.NSites
                        scatter(pos0(ii)+ucrs,zeros(nuc,1),'LineWidth',1.5);
                    end
                    plot(vs,zeros(2,1),'r')
                case 2
                    for ii = 1:uc.NSites
                        posii = num2cell(repmat(pos0(ii,:),nuc,1)+ucrs,1);
                        scatter(posii{:},'LineWidth',1.5);
                    end
                    for ii = 1:2
                        b0 = dec2bin(ii-1,1);
                        ind = @(ib) bin2dec([b0(1:end-ib+1),'0',b0(end-ib+2:end)])+[0,2^(ib-1)]+1;
                        plot(vs(ind(1),1),vs(ind(1),2),'r')
                        plot(vs(ind(2),1),vs(ind(2),2),'b')
                    end
                case 3
                    for ii = 1:uc.NSites
                        posii = num2cell(repmat(pos0(ii,:),nuc,1)+ucrs,1);
                        scatter3(posii{:},'LineWidth',1.5);
                    end
                    for ii = 1:4
                        b0 = dec2bin(ii-1,2);
                        ind = @(ib) bin2dec([b0(1:end-ib+1),'0',b0(end-ib+2:end)])+[0,2^(ib-1)]+1;
                        plot3(vs(ind(1),1),vs(ind(1),2),vs(ind(1),3),'r')
                        plot3(vs(ind(2),1),vs(ind(2),2),vs(ind(2),3),'g')
                        plot3(vs(ind(3),1),vs(ind(3),2),vs(ind(3),3),'b')
                    end
            end
            hold off
            axis tight
            axis equal
            if ~isempty(uc.Labels)
                legend(uc.Labels,'Location','eastoutside')
            end
            set(groot,'defaultAxesColorOrder',c0);
        end
        
%         function tf = isCompatible(uc)
%         end
    end
    
end

