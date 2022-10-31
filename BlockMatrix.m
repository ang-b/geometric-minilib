classdef BlockMatrix < handle
    %BLOCKMATRIX Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = public, SetAccess = private)
        rowBlocks
        columnBlocks
%         rows
%         columns
        rowSizes
        columnSizes
    end
    
    properties (Access = private)
        data; 
    end
    
    methods
        function self = BlockMatrix(varargin)        
            self.rowBlocks = 0;
            self.columnBlocks = 0;
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
                 if isempty(self.data{i,j})
                     if self.rowSizes(i) == 0 && self.columnSizes(j) == 0
                        self.data{i,j} = block;
                        self.rowSizes(i) = blkRows;
                        self.columnSizes(j) = blkCols;
                     else
                         if (self.rowSizes(i) == blkRows && self.columnSizes(j) == blkCols) || ...
                            (self.rowSizes(i) == blkRows && self.columnSizes(j) == 0) || ...
                            (self.rowSizes(i) == 0 && self.columnSizes(j) == blkCols)
                              self.data{i,j} = block;
                         else
                             error('addBlock:cat','Cannot add block of different size');
                         end
                     end
                 end
            end % else it needs to grow
            
            if i > self.rowBlocks 
                self.rowBlocks = i;
                if j > self.columnBlocks
                    self.columnBlocks = j;     
                end        
            end           
        end

        function blk = getBlock(self, i, j)
            if i > length(self.rowSizes) || j > length(self.columnSizes)
                error("out of bounds");
            end
            blk = self.data{i,j};
        end
        
        function M = toMatrix(self)
            for i = 1:self.rowBlocks
                for j = 1:self.columnBlocks
                    if isempty(self.data{i,j})
                        self.data{i,j} = zeros(self.rowSizes(i), self.columnSizes(j));
                    end
                end
            end
            M = cell2mat(self.data);
        end
    end
    
    methods (Access = private)
        
    end
        
end

