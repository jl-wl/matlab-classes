classdef HybridStructure < matlab.mixin.SetGet
    
    properties
        Parts % struct array with fields:
              %   label (string);
              %   tbham (TBHamiltonian);
              %   latps (struct): for values in latps.ns, 0 means infinite,
              %                   -1(-2) means semi-infinite in
              %                   negative(positive) direction (TODO);
        NParts
        Connects % struct array with fields:
                 %   part1 (string): label for part 1;
                 %   ucrs1 (column list of row vectors): unit cell rs at
                 %                               the interface in part 1;
                 %   part2 (string): label for part 2;
                 %   ucrs2 (column list of row vectors): unit cell rs at
                 %                               the interface in part 2;
                 %   tun12 (full tunneling matrix from 2 to 1 at the interface,
                 %        or, function handle generating tunneling matrix)
        DoFs       % cell array of DoF for all parts
        Positions  % cell array of Position for all parts
        Partitions % cell array of Partition for all partsf
        DimHam % dimension of Hamiltonian for all parts
               %   first column for total, second for one unit cell
    end
    
    methods
        function obj = HybridStructure(parts, connects)
        % parts is struct array, or cell array with each row containing:
        %   label, tbham, latps
        % connects is struct array, or cell array with each row containing:
        %   part1, ucrs1, part2, ucrs2, tun12
            if isstruct(parts)
                obj.Parts = parts;
                obj.NParts = length(parts);
            elseif iscell(parts)
                np = size(parts,1);
                for ip = 1:np
                    obj.Parts(ip).label = parts{ip,1};
                    obj.Parts(ip).tbham = parts{ip,2};
                    obj.Parts(ip).latps = parts{ip,3};
                end
                obj.NParts = np;
            else
                error('Invalid input parts!')
            end
            if isstruct(connects)
                obj.Connects = connects;
            elseif iscell(connects)
                for ic = 1:size(connects,1)
                    obj.Connects(ic).part1 = connects{ic,1};
                    obj.Connects(ic).ucrs1 = connects{ic,2};
                    obj.Connects(ic).part2 = connects{ic,3};
                    obj.Connects(ic).ucrs2 = connects{ic,4};
                    obj.Connects(ic).tun12 = connects{ic,5};
                end
            else
                error('Invalid input connects!')
            end
        end
        
        function obj = setParts(obj, varargin)
        % varargin contains the first entry being a number or a label
        %   or a cell array of labels or empty (for all parts),
        %   and the rest being (field, value) pairs
            p = varargin{1};
            if isempty(p)
                ips = 1:length(obj.Parts);
            elseif iscell(p)
                ips = cellfun(@(l)obj.findPart(l),p);
            else
                ips = obj.findPart(p);
            end
            for ii = 1:(length(varargin)-1)/2
                for ip = ips
                    f = varargin{2*ii};
                    v = varargin{2*ii+1};
                    if isa(v,'function_handle')
                        obj.Parts(ip).(f) = v(obj.Parts(ip).(f));
                    else
                        obj.Parts(ip).(f) = v;
                    end
                end
            end
        end
        
        function set.DoFs(obj, dofs)
            obj.DoFs = dofs;
        end
        
        function set.Positions(obj, poss)
            obj.Positions = poss;
        end
        
        function set.Partitions(obj, ptns)
            obj.Partitions = ptns;
        end
        
        function set.DimHam(obj, dh)
            obj.DimHam = dh;
        end
        
        function varargout = eigen(obj, varargin)
            ham_only = true;
            varargout = cell(1,nargout);
            [varargout{:}] = Utility.eigen(obj.fullHamiltonian(ham_only), varargin{:});
        end
        
        function [ham, poss, ptns, dofs] = fullHamiltonian(obj, ham_only)
            if nargout>1 || nargin<2
                ham_only = false;
            end
            ham = [];
            np = obj.NParts;
            dh = zeros(np,2);
            if ~ham_only
                poss = cell(np,1);
                ptns = cell(np,1);
                dofs = cell(np,1);
            end
            for ip = 1:np
                h = obj.Parts(ip).tbham;
                if ~ham_only
                    [hamip, posip, ptnip, dofip] = h.onLattice(obj.Parts(ip).latps);
                    poss{ip} = posip;
                    ptns{ip} = cellfun(@(a)a+sum(dh(1:ip-1,1)),ptnip,'UniformOutput',false);
                    dofs{ip} = dofip;
                else
                    hamip = h.onLattice(obj.Parts(ip).latps);
                end
                ham = blkdiag(ham,hamip);
                dh(ip,:) = [size(hamip,1), h.NDMat];
            end
            if ~ham_only
                obj.set('Positions', poss);
                obj.set('Partitions', ptns);
                obj.set('DoFs', dofs);
            end
            obj.set('DimHam', dh);
            for ic = 1:length(obj.Connects)
                inds1 = obj.findUCIndices(obj.Connects(ic).part1, obj.Connects(ic).ucrs1);
                inds2 = obj.findUCIndices(obj.Connects(ic).part2, obj.Connects(ic).ucrs2);
                V12 = obj.Connects(ic).tun12;
                if isa(V12,'function_handle')
                    ham(inds1,inds2) = V12(obj, ic); %#ok<*AGROW>
                else
                    ham(inds1,inds2) = V12;
                end
                ham(inds2,inds1) = ham(inds1,inds2)';
            end
        end
        
%         function sm = smat(obj, leads) %TODO
%         % To compute scattering matrix with attached leads
%         % leads is struct array with fields:
%             nl = length(leads);
%         end
        
        function ip = findPart(obj, label)
            if isinteger(label)
                ip = label;
            elseif ischar(label)
                ip = find(cellfun(@(l)strcmp(l,label),{obj.Parts.label}));
                if isempty(ip)
                    error(['Part', label, 'is not found!'])
                end
            end
        end
        
        function inds = findUCIndices(obj, label, ucrs)
            if isempty(obj.DoFs)
                ip = obj.findPart(label);
                ns = obj.Parts(ip).latps.ns;
                ns(ns==0) = 1; % treat infinite dimensions as single unit cell
                ucrsip = Utility.ndgridVecs(zeros(size(ns)),ns-1);
                iucs = cellfun(@(v)Utility.findvec(v,ucrsip),mat2cell(ucrs,ones(size(ucrs,1),1))); 
                dH = obj.DimHam(ip,2);
                inds = repelem(iucs,dH,1)*dH+repmat((1-dH:0).',length(iucs),1)+sum(obj.DimHam(1:ip-1,1));
            else
                ip = obj.findPart(label);
                dof = obj.DoFs{ip};
                ucrall = cell2mat(reshape({dof.UnitCell},[],1));
                inds = cellfun(@(v)reshape(Utility.findvec(v,ucrall),[],1),...
                               mat2cell(ucrs,ones(size(ucrs,1),1)),...
                               'UniformOutput',false);
                inds = cell2mat(inds)+sum(obj.DimHam(1:ip-1,1));
            end
        end
    end
            
    methods (Access=private, Hidden=true)
    end
end

