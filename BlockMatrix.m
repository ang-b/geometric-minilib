classdef BlockMatrix < handle
    %BLOCKMATRIX Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = public, SetAccess = private)
        rowSizes
        columnSizes
    end
    
    properties (Access = private)
        data; 
    end
    
    methods
        function self = BlockMatrix(varargin)        
            self.rowSizes = 0;
            self.columnSizes = 0;
            switch (nargin)
                case 0
                    self.data = cell(1);
                case 1 
                    n = varargin{1};
                    self.data = cell(n);
                    self.rowSizes = zeros(n,1);
                    self.columnSizes = zeros(n,1);
                case 2 
                    [m,n] = varargin{:};
                    self.data = cell(m,n);
                    self.rowSizes = zeros(m,1);
                    self.columnSizes = zeros(n,1);
                otherwise
                    error('Wrong input values');
            end
        end
        
        function self = setBlock(self, i, j, block)
            blkRows = size(block, 1);
            blkCols = size(block, 2);
            if i <= size(self.data, 1) && j <= size(self.data, 2)
                if self.blockHasCompatibleRowDimensions(i, block) ...
                        && self.blockHasCompatibleColumnDimensions(j, block)
                    self.data{i,j} = block;
                    self.rowSizes(i) = blkRows;
                    self.columnSizes(j) = blkCols;
                else
                     error('addBlock:cat','Cannot add block of different size');
                end
            end % else it needs to grow
                
        end

        function blk = getBlock(self, i, j)
            if i > length(self.rowSizes) || j > length(self.columnSizes)
                error("out of bounds");
            end
            blk = self.data{i,j};
        end
        
        function s = size(self, varargin)
            if isempty(varargin)
                s = [length(self.rowSizes) length(self.columnSizes)];
            else
                switch varargin{1}
                    case 1
                        s = length(self.rowSizes);
                    case 2
                        s = length(self.columnSizes);
                    otherwise
                        error('BlockMatrix:size', 'Invalid index of dimension');
                end
            end
        end

        function M = toMatrix(self)
            for i = 1:length(self.rowSizes)
                for j = 1:length(self.columnSizes)
                    if isempty(self.data{i,j})
                        self.data{i,j} = zeros(self.rowSizes(i), self.columnSizes(j));
                    end
                end
            end
            M = cell2mat(self.data);
        end
    end
    
    methods (Access = private)
        function flag = blockHasCompatibleColumnDimensions(self, j, block)
            flag = size(block, 2) == self.columnSizes(j) || self.columnSizes(j) == 0;
        end

        function flag = blockHasCompatibleRowDimensions(self, i, block)
            flag = size(block, 1) == self.rowSizes(i) || self.rowSizes(i) == 0;
        end
    end
        
end

