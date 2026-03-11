classdef BinaryTree < handle
    % strict binary trees (2 or 0 descendants)

    properties
        Data
    end

    properties (SetAccess = private)
        Left = []
        Right = []
    end

    methods
        function t = BinaryTree(a, l, r)
            if nargin==0
                a = [];
            end
            t.Data = a;
            if nargin>=3
                t.Left = l;
                t.Right = r;
            end
        end

        function setData(t, data)
            t.Data = data;
        end

        function setLR(t, l, r)
            t.Left = l;
            t.Right = r;
        end

        function tf = isLeaf(t)
            tf = isempty(t.Left);
        end

        function lvs = leaves(t)
            if t.isLeaf
                lvs = {t};
            else
                lvs = [t.Left.leaves, t.Right.leaves];
            end
        end

        function ins = inodes(t) % internal nodes ordered same as "("
            if t.isLeaf
                ins = {};
            else
                ins = [{t}, t.Left.inodes, t.Right.inodes];
            end
        end

        function t = setLeaves(t, ldata)
        % cdata: cell array of length == number of leaves
            lvs = t.leaves;
            arrayfun(@(ii) lvs{ii}.setData(ldata{ii}), 1:length(lvs));
        end

        function dt = upData(t, fn)
        % To update data by applying function fn pairwise, starting from leaves.
        % Note: output of fn must be able to serve as input of fn
            if t.isLeaf
                dt = t.Data;
            else
                dt = fn(t.Left.upData(fn), t.Right.upData(fn));
                t.Data = dt;
            end
        end

        function jt = join(t, t2, fa)
            if nargin>2
                a = fa(t, t2);
            else
                a = [];
            end
            jt = BinaryTree(a, t, t2);
        end

        function [tn, subs] = rot(t, n) % ~ associative move (see assoMove)
            ins = t.inodes;
            if isempty(ins) || n>length(ins) || ins{n}.Left.isLeaf
                tn = [];
                subs = [];
                return
            end
            a = ins{n}.Left.Left;
            b = ins{n}.Left.Right;
            c = ins{n}.Right;
            ancs = BinaryTree.findAncs(t.par, n);
            tn = a.join(b.join(c));
            for ia = size(ancs,2):-1:1
                if ancs(2,ia)<0
                    tn = tn.join(ins{ancs(1,ia)}.Right);
                else
                    tn = ins{ancs(1,ia)}.Left.join(tn);
                end
            end
            subs = {a; b; c};
        end

        function d = depth(t)
            if t.isLeaf
                d = 1;
            else
                d = 1 + max(t.Left.depth, t.Right.depth);
            end
        end

        function w = width(t)
            w = length(t.leaves);
        end

        function s = rep(t, fstr) % text representation
            if nargin<2
                fstr = @string;
            end
            s0 = fstr(t.Data);
            if isempty(s0)
                % warning('Empty data string. Use "X" by default.')
                s0 = "X";
            end
            if t.isLeaf
                s = s0;
            else
                s = pad([string(t.Left.rep(fstr));...
                         string(t.Right.rep(fstr))]);
                s(1) = s0 + " " + s(1);
                s = pad(s, 'left');
            end
        end

        function s = par(t, fstr) % parenthesis representation
            if nargin<2
                fstr = @(a) "";
            end
            if t.isLeaf
                s = fstr(t.Data);
            else
                s = "(" + t.Left.par(fstr) + "," + t.Right.par(fstr) + ")";
            end
        end

        function s = parfull(t, fstr) % par w/ internal data listed behind
            if nargin<2
                fstr = @string;
            end
            s = t.par(fstr);
            s = s + join(cellfun(@(in)fstr(in.Data), t.inodes));
        end
    end

    methods(Static)
        function m = matchPar(par)
        % To match "(", ")", and "," with corresponding indices put in
        % same column, ordered by the occurrence of "(" from left
            p = char(par);
            if ~isequal(count(p,"("), count(p,")"), count(p,","))
                error('Numbers of "(", ")", and "," do not match!')
            end
            inds = [strfind(p,"("); strfind(p,")"); strfind(p,",")];
            n = size(inds,2);
            if n<=1
                m = inds;
                return
            end
            [bs, ix] = sort(reshape(inds(1:2,:),1,[])); % all bracket indices
            ss = reshape([ones(1,n);-ones(1,n)],1,[]); % +1/-1 for (/)
            ss = ss(ix); % ss stores signs corresponding to bs
            m = zeros(3,n);
            for ii = 1:n
                bl = inds(1,ii);
                bsii = bs(bs>=bl);
                ssii = ss(bs>=bl);
                br = bsii(find(~cumsum(ssii),1));
                % cumsum of ss from bl on hits 0 at br
                cm = inds(3,arrayfun(@(ic) ic<br && sum(ssii(bsii<ic))==1, inds(3,:)));
                % there is exactly one more ( than ) from ( to , for THE ,
                if isempty(br) || length(cm)~=1
                    error('Invalid input string format!')
                end
                m(:,ii) = [bl; br; cm];
            end
        end

        function t = genByPar(par, fa)
        % To generate binary tree by parenthesis representation
            par = string(par);
            if nargin<2
                fa = @str2num;
            end
            if ~startsWith(par,"(") % leaf
                t = BinaryTree(fa(par));
                return
            end
            m = BinaryTree.matchPar(par);
            pl = string(extractBetween(par,m(1,1)+1,m(3,1)-1));
            pr = string(extractBetween(par,m(3,1)+1,m(2,1)-1));
            t = BinaryTree([], BinaryTree.genByPar(pl), BinaryTree.genByPar(pr));
        end

        function s = initPar(w) % initial pattern like (((,),),)
            s = "";
            for ii = 2:w
                s = "(" + s + ",)";
            end
        end

        function ss = allPar(w) % all patterns with w leaves
            if w==1
                ss = "";
                return
            end
            ss = cell2mat(reshape(...
                    arrayfun(@(wr) join(Utility.cartprod(...
                    "(",BinaryTree.allPar(w-wr),",",BinaryTree.allPar(wr),")"),""),...
                    1:w-1, 'UniformOutput', false),...
                    [], 1)); % output string array in column
        end

        function s = dataPar(par, data, fstr) % insert data to par
            m = BinaryTree.matchPar(par);
            if length(data)~=size(m,2)+1
                error('Length does not match!')
            end
            if nargin<3
                fstr = @num2str;
            end
            if isnumeric(data)
                data = num2cell(data);
            end
            s = BinaryTree.genByPar(par).setLeaves(data).par(fstr);
        end

        function d = extractData(par, fa) % extract data from par
            par = string(par);
            if nargin<2
                fa = @str2num;
            end
            d = arrayfun(fa, extract(par,regexpPattern('[^(,)]')));
        end

        function as = findAncs(par, n) % ancestors of n-th (
            m = BinaryTree.matchPar(par);
            if n<2 || n>size(m,2)
                as = [];
                return
            end
            ancs = find(m(2,1:n-1)>m(2,n));
            poss = sign(m(3,n)-m(3,ancs)); % -1/+1 for Left/Right
            as = [ancs; poss];
        end

        function [p, subs, sup] = assoMove(par, n)
        % To do associator move ...((_,_),_)... -> ...(_,(_,_))...
        % where the outer ( is the overall n-th (
            par = string(par);
            m = BinaryTree.matchPar(par);
            w = size(m,2)+1;
            if w<3 || n<1 || n>w-2 || m(3,n)<m(1,n+1)
                p = []; % no new tree can be generated
                subs = [];
                sup = [];
                return
            end
            % originally ((a,b),c)
            a = string(extractBetween(par, m(1,n+1)+1, m(3,n+1)-1));
            b = string(extractBetween(par, m(3,n+1)+1, m(2,n+1)-1));
            c = string(extractBetween(par, m(3,n)+1, m(2,n)-1));
            % change to (a,(b,c))
            p = string(replaceBetween(par, m(1,n), m(2,n),...
                "(" + a + ",(" + b + "," + c + "))"));
            subs = [a; b; c];
            sup = string(extractBefore(par,m(1,n)) +...
                    "___" + extractAfter(par,m(2,n)));
        end

        function ps = assoNbs(par) % neighbors by one assoMove
            ps = [];
            for ii = size(BinaryTree.matchPar(par),2):-1:1
                ps = [ps; BinaryTree.assoMove(par,ii)]; %#ok<AGROW>
            end
        end

        function G = associahedron(w, show_plot)
            nodes = BinaryTree.allPar(w);
            nn = length(nodes);
            map = containers.Map(nodes, 1:nn);
            s = [];
            t = [];
            for ii = 1:nn
                nbs = BinaryTree.assoNbs(nodes(ii));
                s = [s; repelem(ii,length(nbs),1)]; %#ok<AGROW>
                t = [t; arrayfun(@(s)map(s),nbs)]; %#ok<AGROW>
            end
            G = digraph(s,t,[],nodes);
            if nargin>1 && show_plot
                figure
                plot(G,'Layout','force3')
            end
        end

        function test
            p = "(1,((2,3),4))";
            t = BinaryTree.genByPar(p);
            disp(p)
            disp(t.rep(@string))
            t.upData(@plus);
            disp("updated tree:")
            disp(t.rep(@string))

            p1 = "((1,((2,3),4)),((5,6),7))";
            t1 = BinaryTree.genByPar(p1);
            t1.upData(@(a,b)string(a)+string(b))
            disp(t1.rep(@string))
            ins = t1.inodes;
            disp('The internal nodes are:')
            celldisp(cellfun(@(c)c.rep(@string), ins, 'UniformOutput', false))
            ii = 4;
            disp(['The ancestors of the ', num2str(ii), '-th internal node:'])
            disp(BinaryTree.findAncs(p1,ii))
            ii = 1;
            [t1r, subs] = t1.rot(ii);
            disp(['Rotation under the ', num2str(ii), '-th internal node:'])
            disp(t1r.rep)
            disp('The three subtrees involved are:')
            celldisp(cellfun(@(c)c.rep(@string), subs, 'UniformOutput', false))

            p2 = '((,),((,),))';
            disp(['The neighbors of ', p2, ' are:'])
            disp(BinaryTree.assoNbs(p2))

            p2d = '((1,2),((3,4),5))';
            ii = 3;
            [pnew, subs, sup] = BinaryTree.assoMove(p2d,ii);
            disp(['The associative move in the ', num2str(ii), '-th',...
                ' "(" of ', p2d, ' results in:'])
            disp(pnew)
            disp('The three subtrees involved are:')
            disp(subs)
            disp('The supertree is:')
            disp(sup)

            w = 5;
            BinaryTree.associahedron(w, true);
            disp(['Showing the associahedron of rank ', num2str(w)])
        end
    end
end