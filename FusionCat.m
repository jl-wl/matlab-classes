classdef FusionCat
    % unitary fusion category

    properties
        SimpleObjs % struct array; first is always tensor unit
        NSimpleObjs
        FusionRules
        Associators
    end

    methods
        function C = FusionCat(sobjs, rules, assos, preset)
            arguments
                sobjs       struct = struct('label','1','PFdim',1)
                rules       double = 1
                assos       cell = {1}
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
            nobj = C.NSimpleObjs;
            assos = cell(nobj,nobj,nobj,nobj);
            iis = Utility.cartprod(1:nobj,1:nobj,1:nobj);
            for ii = num2cell(iis.')
                outs = cellfun(@(t)t.Data, C.fuse(BinaryTree.dataPar('((,),)',ii)));
                dims = arrayfun(@(ii)nnz(outs==ii), 1:nobj);
                for io = 1:nobj
                    if dims(io)==0
                        continue
                    end
                    if isempty(ca)
                        indc = [];
                    else
                        indc = Utility.findvec([ii{:},io], cell2mat(ca(:,1:4)));
                    end
                    if isempty(indc)
                        assos{ii{:},io} = eye(dims(io));
                    else
                        assos{ii{:},io} = ca{indc,5};
                    end
                end
            end
        end

        function ts = fuse(C, par)
        % To generate fusion channels (as BinaryTree) from parenthesis rep
            nobj = C.NSimpleObjs;
            rules = C.FusionRules;
            par = string(par);

            objs = cellfun(@(t)t.Data, BinaryTree.genByPar(par).leaves);
            if any(~ismember(objs, 1:nobj))
                error('Invalid simple object!')
            end

            mbs = BinaryTree.matchPar(par); % matched brackets
            nb = size(mbs,2); % number of brackets
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
            nobj = C.NSimpleObjs;
            rules = C.FusionRules;
            assos = C.Associators;
            par = string(par);

            if nargin<3 || isempty(n) % no move required
                outs = cellfun(@(t)t.Data, C.fuse(par));
                ms = arrayfun(@(ii)eye(nnz(outs==ii)), 1:nobj, 'UniformOutput', false);
                return
            end
            
            [outpar, subs] = BinaryTree.assoMove(par, n);
            if isempty(outpar) % not movable
                ms = {};
                return
            end
            
            ms = cell(nobj,1);
            if n==1
                ams = C.Fmove(subs(1));
                bms = C.Fmove(subs(2));
                cms = C.Fmove(subs(3));
                ias = find(cellfun(@(m)~isempty(m), ams));
                ibs = find(cellfun(@(m)~isempty(m), bms));
                ics = find(cellfun(@(m)~isempty(m), cms));
                iis = Utility.cartprod(ias,ibs,ics,1:nobj);
                for ii = num2cell(iis.')
                    F = assos{ii{:}};
                    if ~isempty(F)
                        ms{ii{4}} = blkdiag(ms{ii{4}},...
                            Utility.kronall(F, ams{ii{1}}, bms{ii{2}}, cms{ii{3}}));
                    end
                end
            else
                m = BinaryTree.matchPar(par);
                lpar = extractBetween(par,m(1,1)+1,m(3,1)-1);
                rpar = extractBetween(par,m(3,1)+1,m(2,1)-1);
                if m(3,n)<m(3,1)
                    nl = n-1;
                    nr = [];
                else
                    nr = n-nnz(m(1,:)<m(3,1));
                    nl = [];
                end
                lms = C.Fmove(lpar, nl);
                rms = C.Fmove(rpar, nr);
                ils = find(cellfun(@(m)~isempty(m), lms));
                irs = find(cellfun(@(m)~isempty(m), rms));
                iis = Utility.cartprod(ils,irs,1:nobj);
                for ii = num2cell(iis.')
                    if rules(ii{:})>0
                        ms{ii{3}} = blkdiag(ms{ii{3}}, kron(lms{ii{1}}, rms{ii{2}}));
                    end
                end
            end
        end

        function tf = checkPantagon(C, tol)
            if nargin<2
                tol = 1e-12;
            end
            nobj = C.NSimpleObjs;
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
                mss = [C.Fmove(dp0,1), C.Fmove(dp0,2), C.Fmove(dp2,1),...
                       C.Fmove(dp3,2), C.Fmove(dp1,1)];
                for io = 1:nobj
                    if isempty(mss{io,1})
                        if any(cellfun(@(m)~isempty(m),mss(io,2:5)))
                            tf = false;
                            disp('Empty matrices do not match!')
                            return
                        end
                    else
                        if norm(mss{io,4}*mss{io,3}*mss{io,2}-mss{io,5}*mss{io,1})>tol
                            tf = false;
                            disp([ii{:}, io])
                            return
                        end
                    end
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