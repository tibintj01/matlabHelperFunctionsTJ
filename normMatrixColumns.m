function nM=normMatrixColumns(M)
    n = sqrt(sum(M.^2,1)); % Compute norms of columns
      nM = bsxfun(@rdivide,M,n); % Normalize M