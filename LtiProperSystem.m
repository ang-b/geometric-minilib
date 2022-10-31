classdef LtiProperSystem < handle
    %LTIPROPERSYSTEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = public, SetAccess = private)
        A
        B
        C
        E
        n
        p
        algebraicProperties
    end

    properties (Access = private)
        outputKernel
        inputSpace
        infUnobserabilitySubspace
        T
        Gamma
        n1
        p1
        G
    end
    
    methods (Access = public)
        function self = LtiProperSystem(A,B,C,E)
            %LTIPROPERSYSTEM Construct an instance of this class
            %   Detailed explanation goes here
            [n,m,p,q,err] = LtiProperSystem.validateConstruction(A,B,C);
            if ~isempty(err), error(err); end
            self.A = A;
            self.B = B;
            self.C = C;
            self.E = E;

            self.n = n;
            self.p = p;
            self.inputSpace = orth(self.B);
            self.outputKernel = null(self.C);

            self.algebraicProperties.outputRedundant = p > m;
            self.algebraicProperties.disturbanceDecomposable = q <= p - m;
            self.algebraicProperties.observable = rank(obsv(A, C)) == n;
            % TODO 
            self.algebraicProperties.leftInvertible = true;
        end
        
        function W = computeInfACinvariantSuspaceContainingB(self)
            W = zeros(self.n);
            for i = 1:self.n
                P = null(W');
                Tw = null([P; self.C]);
                W = [self.B, self.A*Tw];
            end
        end

        function S = computeInfUnobservabilitySubspaceContainingB(self)
            if isempty(self.infUnobserabilitySubspace)
                S = eye(self.n);
                W = self.computeInfACinvariantSuspaceContainingB();
                for i = 1:self.n
                    P = null(S')';
                    Ts = null([P*self.A; self.C]);
                    S = [W, Ts];
                end
                self.infUnobserabilitySubspace = S;
            else
                S = self.infUnobserabilitySubspace;
            end
        end

        function [T, Gamma, n1, p1] = computeUnobSubspaceTransform(self)
            if isempty(self.infUnobserabilitySubspace)
                self.computeInfUnobservabilitySubspaceContainingB();
            end

            if isempty(self.T)
                T2 = self.infUnobserabilitySubspace;
                T1 = null(T2');
                self.n1 = self.n - size(T1, 1);
                self.T = [T1, T2];
                G2 = orth(self.C*T2);
                self.p1 = self.p - size(G2, 2);
                self.Gamma = [null(G2'), G2]; 
            end

            T = self.T;
            Gamma = self.Gamma;
            n1 = self.n1;
            p1 = self.p1;
        end

        function G = designGainLQR(self, Q,R)
            if isempty(self.T)
                self.computeUnobSubspaceTransform();
            end
            
            S = self.infUnobserabilitySubspace;
            % canonical projection onto X/S
            % recall that T1 = null(S');
            P = self.T(:,1:self.n1)'; 

            Ga = -P \ P*self.A*S / (self.C * S);
            Ag = self.A + Ga * self.C;
            Aa = P * Ag / P;
            
            Sc = orth([S self.outputKernel]);
            H = (self.C' \ null(Sc'));
            Ca = H * self.C / P;

            Db = dlqr(Aa.', -Ca.', Q, R).';
            self.G = Ga + P \ Db * H;
            G = self.G;
        end

        function [A,B,C,E,n1,p1] = getUnobSubspaceDecomposition(self)
            if isempty(self.G)
                error("Design a decopuling gain first");
            end

            AG = self.A + self.G * self.C;
            AG = AG .* (abs(AG) > eps);

            A = self.T \ AG * self.T;
            B = self.T \ self.B;
            C = self.Gamma \ self.C * self.T;
            E = self.T \ self.E;

            n1 = self.n1;
            p1 = self.p1;
        end
    end

    methods (Static, Access = private)
        function [n,m,p,q,err] = validateConstruction(A,B,C,E)
            err = [];
            [n, ncols] = size(A);
            if n ~= ncols, err = "A must be a square matrix"; end
            [nb, m] = size(B);
            if nb ~= n, err = "B must have the same number of rows as A"; end
            [ne, q] = size(E);
            if ne ~= n, err = "E must have the same number of rows as A"; end
            [p, nc] = size(C, 2);
            if nc ~= n, err = "C must have the same number of columns as A"; end
            if rank(C) ~= p, err = "Please use full row-rank C"; end
            if rank(B) ~= m, err = "Please use full column-rank B"; end
        end
    end
end

