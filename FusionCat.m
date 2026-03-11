classdef FusionCat
    % unitary fusion category

    properties
        SimpleObjs % struct array; first is always tensor unit
        NSimpleObjs
        FusionRules % (a, b, ab)
        Associators % (a, b, c, abc, bc, ab)
    end

    methods
        function C = FusionCat(sobjs, rules, assos, preset)
            arguments
                sobjs       struct = struct('label','1','PFdim',1)
                rules       double = 1
                assos       double = 1
                preset.name {mustBeMember(preset.name, ...
                            {'tc','fib','is'})}
            end

            C.SimpleObjs = sobjs;
            C.NSimpleObjs = length(C.SimpleObjs);
            C.FusionRules = rules;
            C.Associators = assos;
            
            if isfield(preset,'name')
                switch preset.name
                    case 'tc'
                        C.SimpleObjs = struct('label', {'1';'e';'m';'f'},...
                                              'PFdim', {1;1;1;1});
                        C.NSimpleObjs = 4;
                        C.FusionRules = FusionCat.convertRules(...
                                        {2, 2, 1; 2, 3, 4; 2, 4, 3;...
                                         3, 3, 1; 3, 4, 2;...
                                         4, 4, 1});
                        C.Associators = C.convertAssos({});
                    case 'fib'
                        C.SimpleObjs = struct('label', {'1';'tau'},...
                                              'PFdim', {1;(1+sqrt(sym(5)))/2});
                        C.NSimpleObjs = 2;
                        C.FusionRules = FusionCat.convertRules(...
                                        {2, 2, [1,2]});
                        r = (sqrt(5)-1)/2;
                        C.Associators = C.convertAssos(...
                                        {2, 2, 2, 2, [r,sqrt(r);sqrt(r),-r]});
                    case 'is'
                        C.SimpleObjs = struct('label', {'1';'psi';'sigma'},...
                                              'PFdim', {1;1;sqrt(sym(2))});
                        C.NSimpleObjs = 3;
                        C.FusionRules = FusionCat.convertRules(...
                                        {2, 2, 1; 2, 3, 3; 3, 3, [1,2]});
                        C.Associators = C.convertAssos(...
                                        {2, 3, 2, 3, -1;...
                                         3, 2, 3, 2, -1;...
                                         3, 3, 3, 3, [1,1;1,-1]/sqrt(2)});
                end
            end
        end

        function assos = convertAssos(C, ca)
        % To convert associators from input format (cell array of 5 cols)
        % into storage format (double array of 6 dimensions)
            nobj = C.NSimpleObjs;
            assos = zeros(nobj,nobj,nobj,nobj,nobj,nobj);
            iis = Utility.cartprod(1:nobj,1:nobj,1:nobj);
            for ii = num2cell(iis.')
                p0 = BinaryTree.dataPar('((,),)',ii);
                pn = BinaryTree.assoMove(p0, 1);
                t0s = C.fuse(p0);
                tns = C.fuse(pn);
                out0s = cellfun(@(t)t.Data, t0s);
                outns = cellfun(@(t)t.Data, tns);
                if ~isequal(sort(out0s),sort(outns))
                    error('Inconsistent fusion rules!')
                end
                for io = 1:nobj
                    iabs = sort(cellfun(@(t)t.Left.Data, t0s(out0s==io)));
                    if isempty(iabs) % no fusion channel
                        continue
                    end
                    ibcs = sort(cellfun(@(t)t.Right.Data, tns(outns==io)));
                    if ~isempty(ca)
                        ind = Utility.findvec([ii{:},io], cell2mat(ca(:,1:4)));
                    end
                    if isempty(ca) || isempty(ind) % no nontrivial associator
                        F = eye(length(iabs));
                        % NOTE if length(iabs)>1, this CAN BE WRONG 
                    else
                        F = ca{ind,5};
                    end
                    assos(ii{:},io,ibcs,iabs) = F;
                end
            end
        end

        function ts = fuse(C, par)
        % To generate fusion channels (as BinaryTree) from parenthesis rep
            nobj = C.NSimpleObjs;
            rules = C.FusionRules;
            par = string(par);

            mbs = BinaryTree.matchPar(par); % matched brackets
            nb = size(mbs,2); % number of brackets
            ins = BinaryTree.extractData(par); % input objects
            if length(ins)~=nb+1 || any(~ismember(ins, 1:nobj))
                error('Invalid simple object!')
            end

            if nb<1
                ts = {BinaryTree(str2double(par))};
            elseif nb==1
                l = str2double(extractBetween(par,mbs(1)+1,mbs(3)-1));
                r = str2double(extractBetween(par,mbs(3)+1,mbs(2)-1));
                o = find(rules(l, r, :));
                ts = arrayfun(@(ii) BinaryTree(ii, BinaryTree(l), BinaryTree(r)), o,...
                              'UniformOutput', false);
            else
                tls = C.fuse(extractBetween(par,mbs(1,1)+1,mbs(3,1)-1));
                trs = C.fuse(extractBetween(par,mbs(3,1)+1,mbs(2,1)-1));
                ts = {};
                for il = 1:length(tls)
                    tl = tls{il};
                    for ir = 1:length(trs)
                        tr = trs{ir};
                        o = find(rules(tl.Data, tr.Data, :));
                        ti = arrayfun(@(ii) BinaryTree(ii, tl, tr), o,...
                                      'UniformOutput', false);
                        ts = [ts; ti]; %#ok<AGROW>
                    end
                end
            end
        end

        function ms = Fmove(C, par, n)
        % To generate F-matrices with basis for associative move under
        % n-th internal node (see also BinaryTree.assoMove)
        % ms is (nobj,3) cell array with 1st and 2nd columns for row and col
        % basis as cellarray of BinaryTrees and 3rd column for F-matrix
            nobj = C.NSimpleObjs;
            assos = C.Associators;
            par = string(par);
            
            [pn, subs] = BinaryTree.assoMove(par, n);
            if isempty(pn) % not movable
                ms = {};
                return
            end
            na = count(subs(1),","); % # of internal-nodes in a

            ftdata = @(ts) cellfun(@(t)t.Data, ts);
            % find data of each tree in cell array

            t0s = C.fuse(par);
            out0s = ftdata(t0s);
            tns = C.fuse(pn);
            outns = ftdata(tns);
            if ~isequal(sort(out0s),sort(outns))
                error('Something wrong with channels...')
            end

            ms = cell(nobj,3);
            ms(:,1) = arrayfun(@(ii)tns(outns==ii), (1:nobj)', 'UniformOutput', false); % row basis
            ms(:,2) = arrayfun(@(ii)t0s(out0s==ii), (1:nobj)', 'UniformOutput', false); % col basis
            findata = @(io,rc) cellfun(@(t) ftdata(t.inodes), ms{io,rc}, 'UniformOutput', false);
            % find internal node data of each tree in cell array ms{io,rc}
            del = @(a,n) a((1:length(a))~=n);
            for io = 1:nobj
                if isempty(ms{io,1})
                    continue
                end
                id0s = findata(io,2); % inode data of each basis tree
                idns = findata(io,1); % organized in column cell array
                % inode data in id0s has pattern:
                %  ..., d_abc, d_ab, [na], [nb], [nc], ...
                % inode data in idns has pattern:
                %  ..., d_abc, [na], d_bc, [nb], [nc], ...
                % where na is # of internal nodes in a,
                % [na] is array of length na (same applies to b and c);
                % d_abc appears as the overall n-th internal node.
                % F-matrix is finite only between trees whose inode data
                % differ by d_ab and d_bc only
                dim = length(id0s);
                m = zeros(dim);
                for ii0 = 1:dim
                    [~, subs] = ms{io,2}{ii0}.rot(n);
                    abc = num2cell(ftdata(subs));
                    id0 = id0s{ii0};
                    for iin = 1:dim
                        idn = idns{iin};
                        if isequal(del(idn,n+na+1),del(id0,n+1))
                            m(iin,ii0) = assos(abc{:},id0(n),idn(n+na+1),id0(n+1));
                        end
                    end
                end
                ms{io,3} = m;
            end
        end

        function tf = checkPantagon(C)
            nobj = C.NSimpleObjs;
            mul = @FusionCat.mmul;
            iseq = @FusionCat.ismeq;
            iis = Utility.cartprod(1:nobj,1:nobj,1:nobj,1:nobj);
            p0 = BinaryTree.initPar(4);
            p1 = BinaryTree.assoMove(p0,1);
            p2 = BinaryTree.assoMove(p0,2);
            p3 = BinaryTree.assoMove(p2,1);
            % BinaryTree.assoMove(p3,2) == BinaryTree.assoMove(p1,1)
            for ii = num2cell(iis.')
                dp0 = BinaryTree.dataPar(p0,ii);
                dp1 = BinaryTree.dataPar(p1,ii);
                dp2 = BinaryTree.dataPar(p2,ii);
                dp3 = BinaryTree.dataPar(p3,ii);
                if ~iseq(mul(C.Fmove(dp3,2), mul(C.Fmove(dp2,1), C.Fmove(dp0,2))),...
                         mul(C.Fmove(dp1,1), C.Fmove(dp0,1)))
                    tf = false;
                    return
                end
            end
            tf = true;
        end
    end

    methods(Static)
        function rules = convertRules(ca)
            nobj = max([ca{:,1:2}]);
            rules = zeros(nobj,nobj,nobj);
            rules(1,:,:) = eye(nobj);
            rules(:,1,:) = eye(nobj);
            for ii = ca.'
                rules(ii{:}) = 1;
                if ii{1}~=ii{2}
                    iit = ii([2,1,3]);
                    rules(iit{:}) = 1;
                end
            end
        end

        function m = mmul(m1, m2) % matrix multiplication w/ tree basis
        % m/m1/m2 are cell array with 1st and 2nd columns for row and col
        % basis as cell array of BinaryTrees and 3rd column for matrix
            if ~isequal(size(m1),size(m2))
                error('Size of input does not match!')
            end
            m = cell(size(m1));
            for ii = 1:size(m1,1)
                if isempty(m1{ii,1})
                    m(ii,:) = cell(1,3);
                    continue
                end
                m1c = cell2mat(cellfun(@(t)t.parfull, m1{ii,2}, 'UniformOutput', false));
                m2r = cell2mat(cellfun(@(t)t.parfull, m2{ii,1}, 'UniformOutput', false));
                p = cell2mat(arrayfun(@(s)strcmp(s,m2r'), m1c, 'UniformOutput', false));
                d = length(m1c);
                if ~isequal(1:d,sort((1:d)*p))
                    error('Basis of input does not match!')
                end
                m(ii,:) = {m1{ii,1}, m2{ii,2}, m1{ii,3}*p*m2{ii,3}};
            end
        end

        function tf = ismeq(m1, m2) % isequal bw matrices w/ tree basis
            if ~isequal(size(m1),size(m2))
                disp('Size of input does not match!')
                tf = false;
                return
            end
            for ii = 1:size(m1,1)
                if isempty(m1{ii,1})
                    if ~all(cellfun(@(a)isempty(a), m1(ii,:))) ||...
                       ~all(cellfun(@(a)isempty(a), m2(ii,:)))
                        disp('Empty entries do not match.')
                        tf = false;
                        return
                    else
                        continue
                    end
                end
                m1r = cell2mat(cellfun(@(t)t.parfull, m1{ii,1}, 'UniformOutput', false));
                m1c = cell2mat(cellfun(@(t)t.parfull, m1{ii,2}, 'UniformOutput', false));
                m2r = cell2mat(cellfun(@(t)t.parfull, m2{ii,1}, 'UniformOutput', false));
                m2c = cell2mat(cellfun(@(t)t.parfull, m2{ii,2}, 'UniformOutput', false));
                pr = cell2mat(arrayfun(@(s)strcmp(s,m2r'), m1r, 'UniformOutput', false));
                pc = cell2mat(arrayfun(@(s)strcmp(s,m2c), m1c', 'UniformOutput', false));
                dr = length(m1r);
                if ~isequal(1:dr,sort((1:dr)*pr))
                    disp('Row basis does not match.')
                    tf = false;
                    return
                end
                dc = length(m1c);
                if ~isequal(1:dc,sort((1:dc)*pc))
                    disp('Column basis does not match.')
                    tf = false;
                    return
                end
                if ~isapprox(m1{ii,3}, pr*m2{ii,3}*pc, "tight")
                    disp('Matrix does not match.')
                    tf = false;
                    return
                end
            end
            tf = true;
        end

        function test
            TC = FusionCat(name='tc');
            disp(TC)
            ttc = TC.fuse('(((1,2),3),4)');
            celldisp(cellfun(@(c)c.rep, ttc, 'UniformOutput', false))

            Is = FusionCat(name='is');
            disp(Is)
            tis = Is.fuse('((3,3),(3,3))');
            celldisp(cellfun(@(c)c.rep, tis, 'UniformOutput', false))

            Fib = FusionCat(name='fib');
            disp(Fib)
            tfib = Fib.fuse('(((2,2),(2,2)),2)');
            celldisp(cellfun(@(c)c.rep, tfib, 'UniformOutput', false))
        end
    end
end