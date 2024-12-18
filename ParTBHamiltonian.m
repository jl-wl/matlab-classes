classdef ParTBHamiltonian < matlab.mixin.SetGet
    % Parametrized Tight-Binding Hamiltonian
    
    properties
        Hoppings % Hoppings is a nhops*3 cell array w/ each row containing
                 %   f : empty, or anonymous function handle
                 %   R : hopping (row) vector
                 %   HR: hopping matrix
        DoF      % degrees of freedom (struct array)
        ParList  % cell array: list of parameters defined in Hoppinggs
        Parameters % struct array: values of parameters
        NDims    % spacial dimension
        NDMat    % number of dimensions of HR
    end
    
    methods
        function h = ParTBHamiltonian(hops, dof)
        % hops is character array (commonly used model) or cell array w/
        %   each row of hops contains f, Rs, HR
        %   f can be empty, par name(str), or anonymous function handle
        %   Rs is a nR*dim array
        %   HR is a single matrix or cell arry of length nR
        %   [NOTE] hermitian conjugate terms will be automatically added if R~=0
        % dof can be struct array or input format of TBHamiltonian.convertDoF
            if ischar(hops)
                if nargin<2
                    [hops, dof] = ParTBHamiltonian.findModel(hops);
                else % the input dof will override the original one
                    hops = ParTBHamiltonian.findModel(hops);
                end
                if isempty(hops)
                    error('Unrecognized model name!')
                end
            elseif nargin<2
                dof = [];
            end
            h.Hoppings = {};
            h.ParList = {};
            h.addHops(hops);
            h.NDims = size(h.Hoppings{1,2},2);
            h.NDMat = size(h.Hoppings{1,3},1);
            if ~isempty(dof)
                if iscell(dof)
                    h.DoF = TBHamiltonian.convertDoF(dof);
                elseif isstruct(dof)
                    h.DoF = dof;
                else
                    error('Invalid dof input!')
                end
            else
                h.DoF = [];
            end
        end
        
        function h = addHops(h, hops)
            for ii = 1:size(hops,1)
                f = hops{ii,1};
                Rs = hops{ii,2};
                HR = hops{ii,3};
                if ischar(f)
                    p = {f};
                    f = str2func(['@(',f,')',f]);
                elseif ~isempty(f) % assume to be anonymous function handle
                    fstr = func2str(f);
                    p = split(strtok(fstr(3:end),')'),',');
                else % f is empty
                    p = {};
                end
                nR = size(Rs,1);
                if iscell(HR) && nR~=length(HR)
                    error('Numbers of hopping vectors and matrices not compatible!')
                end
                for iR = 1:nR
                    if iscell(HR)
                        h.Hoppings = cat(1,h.Hoppings,{f,Rs(iR,:),HR{iR}});
                    else
                        h.Hoppings = cat(1,h.Hoppings,{f,Rs(iR,:),HR});
                    end
                    h.ParList = cat(1,h.ParList,p);
                end
                h.ParList = unique(h.ParList);
            end
        end
        
        function set.Parameters(h, p)
            if isstruct(p)
                h.Parameters = p;
            else
                for ii = 1:length(p)/2
                    h.Parameters.(p{2*ii-1}) = p{2*ii};
                end
            end
        end
        
        function h = setPar(h, varargin)
        % convenient interface to set parameters
            if isscalar(varargin)
                h.set('Parameters',varargin{1});
            else
                h.set('Parameters',varargin);
            end
        end
        
        function hn = double(h0, dofn, transform, add_hops)
        % double number of DoF, e.g. to add SOC or pairing
        % dofn can be full struct array for the new Hamiltonian or cell
        %   array containing only the new dof, e.g. {'Spin',{'u','d'}} or
        %   {'Nambu',{'e','h'}}
        % transform is a function handle taking (f, R, HR, DoF) as input
            hops0 = h0.Hoppings;
            dof0 = h0.DoF;
            n = size(hops0,1);
            dH0 = h0.NDMat;
            hopsn = cell(2*n,3);
            for ii = 1:n
                hopsn(ii,:) = hops0(ii,:);
                hopsn{ii,3} = blkdiag(hopsn{ii,3},zeros(dH0));
                [hopsn{ii+n,:}] = transform(hops0{ii,:},dof0);
                hopsn{ii+n,3} = blkdiag(zeros(dH0),hopsn{ii+n,3});
            end
            
            if iscell(dofn)
                name = dofn{1};
                vals = dofn{2};
                dofn = repmat(dof0,2,1);
                for ii = 1:dH0
                    dofn(ii).(name) = vals{1};
                    dofn(ii+dH0).(name) = vals{2};
                end
            end
            
            hn = ParTBHamiltonian(hopsn,dofn);
            hn.setPar(h0.Parameters);
            if nargin>3
                hn.addHops(add_hops);
            end
        end
        
        function hh = build(h, p, varargin)
        % to build a specific TBHamiltonian from parametrized one
            if nargin>=2
                h.set('Parameters',p);
            end
            pars = h.Parameters;
            hops = h.Hoppings;
            Rs = [];
            HR = [];
            for ii = 1:size(hops,1)
                f = hops{ii,1};
                if isempty(f)
                    v = 1;
                else
                    fstr = func2str(f);
                    p = split(strtok(fstr(3:end),')'),',');
                    if ~iscell(p)
                        p = cellstr(p);
                    end
                    in = cellfun(@(s)pars.(s),p,'UniformOutput',false);
                    v = f(in{:});
                end
                Rii = hops{ii,2};
                HRi = hops{ii,3}*v;
                Rs = cat(1,Rs,Rii);
                HR = cat(3,HR,full(HRi));
                if sum(abs(Rii))>0 || ~ishermitian(HRi) % add hermitian conjugate
                    Rs = cat(1,Rs,-Rii);
                    HR = cat(3,HR,full(HRi'));
                end
            end
            hh = TBHamiltonian(Rs, HR, h.DoF, varargin{:});
        end
    end
    
    methods (Static)
        function [hops, dof] = findModel(name)
            util = Utility();
            [s1, s2, s3, s0] = util.PauliMats{:};
            switch name
                case 'spinless1d' % 1D spinless
                    dof = [];
                    hops = {@(t,mu)2*t-mu, 0, 1; ...
                            @(t)-t,        1, 1};
                case 'spinless2d' % 2D spinless
                    dof = [];
                    hops = {@(t,mu)4*t-mu, [0,0],  1; ...
                            @(t)-t,        eye(2), 1};
                case 'spinless3d' % 3D spinless
                    dof = [];
                    hops = {@(t,mu)6*t-mu, [0,0,0], 1; ...
                            @(t)-t,        eye(3),  1};
                case 'spin1d' % 1D spinful
                    dof = {...
                           {'',1:2,...
                            'Spin',{'u','d'}}};
                    hops = {@(t,mu)2*t-mu, 0, s0; ...
                            @(t)-t,        1, s0; ...
                            @(M,th_M,ph_M)M*sin(th_M)*cos(ph_M), 0, s1; ...
                            @(M,th_M,ph_M)M*sin(th_M)*sin(ph_M), 0, s2; ...
                            @(M,th_M)     M*cos(th_M),           0, s3; ...
                            @(t_SO,th_SO,ph_SO)1i*t_SO*sin(th_SO)*cos(ph_SO), 1, s1; ...
                            @(t_SO,th_SO,ph_SO)1i*t_SO*sin(th_SO)*sin(ph_SO), 1, s2; ...
                            @(t_SO,th_SO)      1i*t_SO*cos(th_SO),            1, s3};
                case 'rashba2d' % 2D with Rashba SOC and Zeeman
                    dof = {...
                           {'',1:2,...
                            'Spin',{'u','d'}}};
                    hops = {@(t,mu)4*t-mu, [0,0],  s0; ...
                            @(t)-t,        eye(2), s0; ...
                            @(M,th_M,ph_M)M*sin(th_M)*cos(ph_M), [0,0], s1; ...
                            @(M,th_M,ph_M)M*sin(th_M)*sin(ph_M), [0,0], s2; ...
                            @(M,th_M)     M*cos(th_M),           [0,0], s3; ...
                            @(t_R)1i*t_R,  [1,0], s2; ...
                            @(t_R)-1i*t_R, [0,1], s1};
                case {'qah2d', 'ci2d'} % quantum anomalous Hall / Chern insulator in 2D
                    dof = [];
                    hops = {@(m,b)4*b-m, [0,0],  s3; ...
                            @(b)-b,      eye(2), s3;...
                            @(a)0.5i*a,  eye(2), {s1,s2}};
                case 'graphene'
                    dof = {...
                           {'A',1};...
                           {'B',2}};
                    hops = {'t', [0,0],  s1; ...
                            't', [1,0],  [0,1;0,0]; ...
                            't', [0,1],  [0,1;0,0]};
                case 'ssc1d' % s-wave superconductor in 1D (ignore spin)
                    dof = {...
                           {'',1:2,...
                            'Nambu',{'e','h'}}};
                    hops = {@(t,mu)2*t-mu, 0, s3; ...
                            @(t)-t,        1, s3;...
                            @(D,ph_D)D*cos(ph_D), 0, s1;...
                            @(D,ph_D)D*sin(ph_D), 0, s2};
                case 'psc1d' % chiral p-wave superconductor in 1D
                    dof = {...
                           {'',1:2,...
                            'Nambu',{'e','h'}}};
                    hops = {@(t,mu)2*t-mu, 0, s3; ...
                            @(t)-t,        1, s3;...
                            @(D,ph_D)0.5i*D*cos(ph_D), 1, s1;...
                            @(D,ph_D)0.5i*D*sin(ph_D), 1, s2};
                case 'dsc1d' % d-wave superconductor in 1D
                    dof = {...
                           {'',1:2,...
                            'Nambu',{'e','h'}}};
                    hops = {@(t,mu)2*t-mu, 0, s3; ...
                            @(t)-t,        1, s3;...
                            @(D,ph_D)-D*cos(ph_D), 1, s1;...
                            @(D,ph_D)-D*sin(ph_D), 1, s2};
                case 'ssc2d' % s-wave superconductor in 2D (ignore spin)
                    dof = {...
                           {'',1:2,...
                            'Nambu',{'e','h'}}};
                    hops = {@(t,mu)4*t-mu, [0,0], s3; ...
                            @(t)-t,        eye(2), s3;...
                            @(D,ph_D)D*cos(ph_D), [0,0], s1;...
                            @(D,ph_D)D*sin(ph_D), [0,0], s2};
                case 'psc2d' % chiral p-wave superconductor in 2D
                    dof = {...
                           {'',1:2,...
                            'Nambu',{'e','h'}}};
                    hops = {@(t,mu)4*t-mu, [0,0],  s3; ...
                            @(t)-t,        eye(2), s3;...
                            @(D,ph_D)0.5i*D*cos(ph_D), eye(2), {s1,s2};...
                            @(D,ph_D)0.5i*D*sin(ph_D), eye(2), {s2,-s1}};
                case 'psc2d_aniso' % anisotropic chiral p-wave sc in 2D
                    dof = {...
                           {'',1:2,...
                            'Nambu',{'e','h'}}};
                    hops = {@(tx,ty,mu)2*tx+2*ty-mu, [0,0], s3; ...
                            @(tx)-tx,                [1,0], s3;...
                            @(ty)-ty,                [0,1], s3;...
                            @(Dx,ph_D)0.5i*Dx*cos(ph_D), [1,0], s1;...
                            @(Dx,ph_D)0.5i*Dx*sin(ph_D), [1,0], s2;...
                            @(Dy,ph_D)0.5i*Dy*cos(ph_D), [0,1], s2;...
                            @(Dy,ph_D)0.5i*Dy*sin(ph_D), [0,1], -s1};
                case 'dxysc' % d(2xy)-wave superconductor in 2D
                    dof = {...
                           {'',1:2,...
                            'Nambu',{'e','h'}}};
                    hops = {@(t,mu)4*t-mu, [0,0],  s3; ...
                            @(t)-t,        eye(2), s3;...
                            @(D,ph_D)(D/2)*cos(ph_D), [1,-1;1,1], {s1,-s1};...
                            @(D,ph_D)(D/2)*sin(ph_D), [1,-1;1,1], {s2,-s2}};
                case 'dx2y2sc' % d(x^2-y^2)-wave superconductor in 2D
                    dof = {...
                           {'',1:2,...
                            'Nambu',{'e','h'}}};
                    hops = {@(t,mu)4*t-mu, [0,0],  s3; ...
                            @(t)-t,        eye(2), s3;...
                            @(D,ph_D)-D*cos(ph_D), eye(2), {s1,-s1};...
                            @(D,ph_D)-D*sin(ph_D), eye(2), {s2,-s2}};
                case 'ti2d'
                    dof = {...
                           {'',1:4,...
                            'Spin',{'u','d'},...
                            'Orbital',{'s','p'}}};
                    hops = {@(mu)-mu,    [0,0],  eye(4);...
                            @(m,b)4*b-m, [0,0],  kron(s0,s3);...
                            'Dsox',      [0,0],  kron(s1,s2);...
                            'Dsoy',      [0,0],  kron(s2,s2);...
                            @(b)-b,      eye(2), kron(s0,s3);...
                            @(a)0.5i*a,  eye(2), {kron(s3,s1),kron(s0,s2)}};
%                             @(a)0.5i*a,  eye(2), {kron(s0,s2),kron(s3,s1)}};
                case 'ti3d'
                    dof = [];
                    beta = kron(s0,s3);
                    alphas = {kron(s1,s1), kron(s2,s1), kron(s3,s1)};
                    hops = {@(mu)-mu,    [0,0,0], eye(4);...
                            @(m,b)6*b-m, [0,0,0], beta;...
                            @(b)-b,      eye(3),  beta;...
                            @(a)0.5i*a,  eye(3),  alphas};
                case 'weyl3d'
                    dof = [];
                    hops = {@(r,m,d)4*d-r*m, [0,0,0],        s0;...
                            @(d)-d,          [1,0,0; 0,1,0], s0;...
                            @(r,f)r*f,       [0,0,1],        s0;...
                            @(a)1i*a,        [1,0,0; 0,1,0], {s1,s2};...
                            @(m,b)4*b-m,     [0,0,0],        s3;...
                            @(b)-b,          [1,0,0; 0,1,0], s3;...
                            'f',             [0,0,1],        s3};
                case 'hoti2d'
                    G = {kron(s3,s1), kron(s3,s2), kron(s3,s3),...
                         kron(s1,s0), kron(s2,s0)};
                    dof = [];
                    hops = {@(m,b)4*b-m, [0,0],  G{3};...
                            @(b)-b,      eye(2), G{3};...
                            @(a)0.5i*a,  eye(2), G(1:2);...
                            @(d) d,      [1,0], G{4};...
                            @(d)-d,      [0,1], G{4}};
                case 'triangle' % spinless 2D
                    dof = [];
                    hops = {'t1',        [1, 0], 1; ...
                            't2',        [-1,1], 1; ...
                            't3',        [0,-1], 1};
                otherwise
                    try
                        f = str2func(name);
                        [hops, dof] = f();
                    catch
                        warning(['Model ',name,' not found!']);
                        dof = [];
                        hops = {};
                    end
            end
        end
        
        function [f, R, HR] = TR1(f0, R0, HR0, ~)
        % time reversal transformation for one spin species
        %   simply complex conjugate is applied
            [tok, rem] = strtok(func2str(f0),')');
            tok = [tok,')'];
            rem = ['conj(',rem(2:end),')'];
            f = str2func([tok,rem]);
            R = R0;
            HR = conj(HR0);
        end
        
        function [f, R, HR] = minusTR2(f0, R0, HR0, dof)
        % time reversal transformation for two spin species, multiplied by
        % a minus sign at last, as used in BdG Hamiltonian
        %   T = i*sigma_y*K is used for the TR operator
            [tok, rem] = strtok(func2str(f0),')');
            tok = [tok,')'];
            rem = ['conj(',rem(2:end),')'];
            f = str2func([tok,rem]);
            R = R0;
            if nargin<4 || isempty(dof)
                UT = kron([0,1;-1,0],eye(size(HR0,1)/2));
            else
                UT = TBHamiltonian.fullOperator(dof, {'Spin',{'u','d'},[0,1;-1,0]});
            end
            HR = UT*conj(HR0)*UT;
        end
    end
    
end

