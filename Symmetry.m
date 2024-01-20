classdef Symmetry
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        function [Js, bas] = am(j, basis)
        % angular momentum operators Js = {Jx, Jy, Jz}
            if nargin<2
                basis = 'spherical';
            end
            twoj = 2*j;
            assert(~mod(twoj,1), 'j must be a nonnegative integer or half-integer')
            d = round(twoj+1);
            m = (twoj:-2:-twoj)'/2;
            Js{3} = Utility.chop(spdiags(m,0,d,d));
            Jplus = spdiags(sqrt((j+m+1).*(j-m)),1,d,d);
            Js{1} = Utility.chop(0.5*(Jplus+Jplus'));
            Js{2} = Utility.chop(-0.5i*(Jplus-Jplus'));
            if ~mod(j,1) && strcmp(basis,'cubic')
                [U, bas] = Symmetry.sph2cub(j);
                for ii=1:3
                    Js{ii} = Utility.chop(U*Js{ii}*U');
                end
                bas = bas(1,:);
            else
                bas = arrayfun(@(mz)['mz=',num2str(mz)],m','UniformOutput',false);
            end
        end
        
        function [U, basis] = sph2cub(L)
        % Return transformation matrix from spherical to cubic harmonics.
        % Notes
        % -----
        % The spherical harmonics basis is assumed to be ordered from +L to -L;
        % the cubic harmics basis is defined from the spherical basis as
        % 
        % Y_L^{ym} = (Y_L^m - (-1)**m Y_L^{-m}) * 1j/sqrt(2), m = L, L-1, ..., 1,
        % Y_L^{z} = Y_L^0,
        % Y_L^{xm} = (Y_L^m + (-1)**m Y_L^{-m}) * 1/sqrt(2), m = 1, 2, ..., L.
        % 
        % This is not exactly the same as in the convention with the spherical
        % basis ordered from -L to +L, e.g. in
        % https://en.wikipedia.org/wiki/Spherical_harmonics#Real_form
            assert((L>=0)&&~mod(L,1), 'L must be a non-negative integer')
            switch L
                case 0
                    basis = {'s'; 'mz=0'};
                case 1
                    basis = {'py', 'pz', 'px'; 'mz=1', 'mz=0', 'mz=-1'};
                case 2
                    basis = {'dxy', 'dyz', 'dz^2', 'dzx', 'dx^2-y^2';...
                             'mz=2', 'mz=1', 'mz=0', 'mz=-1', 'mz=-2'};
                case 3
                    basis = {'fy(3x^2-y^2)', 'fxyz', 'fyz^2', 'fz^3',...
                                   'fxz^2', 'fz(x^2-y^2)', 'fx(x^2-3y^2)';
                             'mz=3', 'mz=2', 'mz=1', 'mz=0', 'mz=-1', 'mz=-2', 'mz=-3'};
                otherwise
                    basis = cat(1,...
                        cat(2,arrayfun(@(n)['y',num2str(n)],L:-1:1,'UniformOutput',false),...
                        'z',arrayfun(@(n)['x',num2str(n)],1:L,'UniformOutput',false)),...
                        arrayfun(@(m)['mz=',num2str(m)],L:-1:-L,'UniformOutput',false));
            end
            if L==0
                U = 1;
            else
                dp = sqrt(0.5)*eye(L);
                dn = sqrt(0.5)*diag((-1).^(1:L));
                U = [1j*dp, zeros(L,1), -1j*flipud(dn);...
                     zeros(1,L), 1, zeros(1,L);...
                     flipud(dp), zeros(L,1), dn];
            end
        end

        function Ys = sphharm(n, phi, theta)
        % computes spherical harmonics of degree n and order m = 0,1,...,n
        % phi : Azimuthal (longitudinal) coordinate; must be in [0, 2*pi]
        % theta : Polar (colatitudinal) coordinate; must be in [0, pi]
            Ys = (-exp(1i*phi)).^((0:n).').*legendre(n,cos(theta),'norm')/sqrt(2*pi);
        end
        
        function w = wigner3j( j123, m123 )
        % Compute the Wigner 3j symbol using the Racah formula. 
        %
        % W = Wigner3j( J123, M123 ) 
        %
        % J123 = [J1, J2, J3].
        % M123 = [M1, M2, M3].
        % All Ji's and Mi's have to be integeres or half integers (correspondingly).
        %
        % According to seletion rules, W = 0 unless:
        %   |Ji - Jj| <= Jk <= (Ji + Jj)    (i,j,k are permutations of 1,2,3)
        %   |Mi| <= Ji    (i = 1,2,3)
        %    M1 + M2 + M3 = 0
        % 
        % Reference: 
        % Wigner 3j-Symbol entry of Eric Weinstein's Mathworld:
        % http://mathworld.wolfram.com/Wigner3j-Symbol.html
        %
        % Inspired by Wigner3j.m by David Terr, Raytheon, 6-17-04
        %  (available at www.mathworks.com/matlabcentral/fileexchange).
        %
        % By Kobi Kraus, Technion, 25-6-08.
        % Updated 1-8-13.
            j1 = j123(1); j2 = j123(2); j3 = j123(3);
            m1 = m123(1); m2 = m123(2); m3 = m123(3);
            % Input error checking
            if any( j123 < 0 )
                error( 'The j must be non-negative' )
            elseif any( rem( [j123, m123], 0.5 ) )
                error( 'All arguments must be integers or half-integers' )
            elseif any( rem( (j123 - m123), 1 ) )
                error( 'j123 and m123 do not match' );
            end
            % Selection rules
            if ( j3 > (j1 + j2) ) || ( j3 < abs(j1 - j2) ) ... % j3 out of interval
               || ( m1 + m2 + m3 ~= 0 ) ... % non-conserving angular momentum
               || any( abs( m123 ) > j123 ) % m is larger than j
                w = 0;
                return
            end
            % Simple common case
            if ~any( m123 ) && rem( sum( j123 ), 2 ) % m1 = m2 = m3 = 0 & j1 + j2 + j3 is odd
                w = 0;
                return
            end
            % Evaluation
            t1 = j2 - m1 - j3;
            t2 = j1 + m2 - j3;
            t3 = j1 + j2 - j3;
            t4 = j1 - m1;
            t5 = j2 + m2;
            tmin = max( 0,  max( t1, t2 ) );
            tmax = min( t3, min( t4, t5 ) );
            t = tmin : tmax;
            w = sum( (-1).^t .* exp( -ones(1,6) * gammaln( [t; t-t1; t-t2; t3-t; t4-t; t5-t] +1 ) + ...
                                     gammaln( [j1+j2+j3+1, j1+j2-j3, j1-j2+j3, -j1+j2+j3, j1+m1, j1-m1, j2+m2, j2-m2, j3+m3, j3-m3] +1 ) ...
                                     * [-1; ones(9,1)] * 0.5 ) ) * (-1)^( j1-j2-m3 );
            % Warnings
            if isnan( w )
                warning( 'MATLAB:Wigner3j:NaN', 'Wigner3J is NaN!' )
            elseif isinf( w )
                warning( 'MATLAB:Wigner3j:Inf', 'Wigner3J is Inf!' )
            end
        end

        function inv_op = inversionPosition(pos)
            pos = pos-repelem((min(pos)+max(pos))/2,size(pos,1),1);
            try
                ind = cellfun(@(p)Utility.findvec(-p,pos),num2cell(pos,2));
            catch
                error('Positions do not satisfy inversion symmetry!')
            end
            inv_op = sparse(1:length(ind),ind,ones(1,length(ind)));
        end
    end
    
end