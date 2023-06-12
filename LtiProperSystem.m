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
        G
        algebraicProperties
        outputKernel
        inputSpace
        inputMap
        infUnobserabilitySubspace
        T
        Gamma
        n1
        p1
        transformedSystem
    end
    
    methods (Access = public)
        function self = LtiProperSystem(A,B,C,E)
            %LTIPROPERSYSTEM Construct an instance of this class
            %   Detailed explanation goes here
            [n,m,p,q,err] = LtiProperSystem.validateConstruction(A,B,C,E);
            if ~isempty(err), error(err); end
            self.A = A;
            self.B = B;
            self.C = C;
            self.E = E;

            self.n = n;
            self.p = p;
            self.inputMap = [];
            self.inputSpace = [];
            self.outputKernel = null(self.C);

            self.algebraicProperties.outputRedundant = p > m;
            self.algebraicProperties.disturbanceDecomposable = q <= p - m;
            self.algebraicProperties.observable = rank(obsv(A, C)) == n;
            % TODO 
            self.algebraicProperties.leftInvertible = true;
        end

        function self = setInputMap(self, imap)
            switch imap
                case {'b', 'B'}
                    self.inputMap = self.B;
                case {'e', 'E'}
                    self.inputMap = self.E;
                otherwise, error("Choice must be either matrix E or B");
            end
            self.inputSpace = orth(self.inputMap);
            self.infUnobserabilitySubspace = [];
            self.T = [];
            self.G = [];
        end

        function W = getInfACinvariantSuspaceContainingInputSpace(self)
            if isempty(self.inputMap) || isempty(self.inputSpace)
                error("Select the input map first");
            end
            W = zeros(self.n);
            for i = 1:self.n
                P = null(W')';
                Tw = null([P; self.C]);
                W = [self.inputMap, self.A*Tw];
            end
        end

        function S = getInfUnobservabilitySubspaceContainingInput(self)
            if isempty(self.infUnobserabilitySubspace)
                S = eye(self.n);
                W = self.getInfACinvariantSuspaceContainingInputSpace();
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

        function [T, Gamma, n1, p1] = getUnobSubspaceTransform(self)
            if isempty(self.infUnobserabilitySubspace)
                self.getInfUnobservabilitySubspaceContainingInput();
            end

            if isempty(self.T)
                T2 = self.infUnobserabilitySubspace;
                T1 = null(T2');
                self.n1 = self.n - size(T2, 2);
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

        function G = placePoles(self, usPoles, usModXPoles)
            if isempty(self.T)
                self.getUnobSubspaceTransform();
            end
            S = self.infUnobserabilitySubspace;
            % canonical projection onto X/S
            % recall that T1 = null(S');
            P = self.T(:,1:self.n1)'; 

            L0 = -P \ P*self.A*S / (self.C * S);
            AL = self.A + L0 * self.C;
            A_S = S \ AL * S;
            C_S = self.C*S;

            L_S = -place(A_S',C_S',usPoles)';
            L0 = L0 + S*L_S;
            
            A_XS = P * AL / P;
            Sc = orth([S null(self.C,'r')]);
            H = (self.C' \ null(Sc'))';
            C_XS = H * self.C / P;
        
            L_XS = -place(A_XS', C_XS', usModXPoles)';
            G = L0 + P \ L_XS * H;
        end

%         function G = designGainLQR(self,Q,R)
%             if isempty(self.T)
%                 self.getUnobSubspaceTransform();
%             end
%             
%             % args Q and R are relative to coordinates x,y
%             Q = self.T'*Q*self.T;
% %             R = self.Gamma' \ R / self.Gamma;
% 
%             S = self.infUnobserabilitySubspace;
%             % canonical projection onto X/S
%             % recall that T1 = null(S');
%             P = self.T(:,1:self.n1)'; 
%             Sc = orth([S self.outputKernel]);
%             
%             H_XS = (self.C' \ null(Sc', 'r'))';
%             H_S = null(H_XS, 'r');
% 
%             G0 = -P \ P*self.A*S / (self.C * S);
%             Ag = self.A + G0 * self.C;
%             
%             % projection on the unobservability subspace
%             A_S = S \ Ag * S;
%             C_S = self.C*S;
%             Q_S = S \ Q * S;
%             R_S = H_S \ R;
% 
%             % TODO: fix this
% %             L_S = dlqr(A_S.', -C_S.', Q_S, R_S).';
% %             G0 = G0 + S*L_S;
% 
%             A_XS = P * self.A / P;
%             Q_XS = P * Q / P;
%                         
%             C_XS = H_XS * self.C / P;
%             % TODO: fix dimensions of this
%             R_XS = H_XS * R;
% 
%             L_XS = dlqr(A_XS.', -C_XS.', Q_XS, R_XS).';
%             self.G = G0 + P \ L_XS * H_XS;
%             G = self.G;
%         end

        function setObserverFeedback(self, G)
            self.G = G;
        end

        function [A,B,C,E,n1,p1] = getUnobSubspaceDecomposition(self)
            if isempty(self.G)
                error("Design a decopuling gain first");
            end
            if isempty(self.T)
                self.getUnobSubspaceTransform();
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

        function sys = getBlockUnobSubspaceDecomposition(self)
            [AT,BT,CT,ET,~,~] = self.getUnobSubspaceDecomposition();

            ABlock = BlockMatrix(2,2);
            ABlock.setBlock(1, 1, AT(1:self.n1, 1:self.n1));
            ABlock.setBlock(2, 1, AT(self.n1+1:self.n, 1:self.n1));
            ABlock.setBlock(2, 2, AT(self.n1+1:self.n, self.n1+1:self.n));
            self.transformedSystem.A = ABlock;
            
            BBlock = BlockMatrix(2,1);
            BBlock.setBlock(2,1, BT(self.n1+1:self.n, :));
            self.transformedSystem.B = BBlock;

            CBlock = BlockMatrix(2,2);
            CBlock.setBlock(1, 1, CT(1:self.p1, 1:self.n1));
            CBlock.setBlock(2, 1, CT(self.p1+1:self.p, 1:self.n1));
            CBlock.setBlock(2, 2, CT(self.p1+1:self.p, self.n1+1:self.n));
            self.transformedSystem.C = CBlock;

            EBlock = BlockMatrix(2,1);
            EBlock.setBlock(1, 1, ET(1:self.n1, :));
            EBlock.setBlock(2, 1, ET(self.n1+1:self.n, :));
            self.transformedSystem.E = EBlock;

            sys = self.transformedSystem;
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
            [p, nc] = size(C);
            if nc ~= n, err = "C must have the same number of columns as A"; end
            if rank(C) ~= p, err = "Please use full row-rank C"; end
            if rank(B) ~= m, err = "Please use full column-rank B"; end
        end
    end
end

