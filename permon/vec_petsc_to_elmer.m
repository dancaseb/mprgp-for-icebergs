args = argv();
if numel(args) >= 1
  inputFile = args{1};
else
  inputFile = 'x_output.bin';
end
if numel(args) >= 2
  outputFile = args{2};
else
  outputFile = 'x_output.dat';
end

b = [];   % vector that will hold data
success = false;

% read
if exist('PetscBinaryRead', 'file') || exist('PetscBinaryRead', 'builtin')
  try
    fprintf('Attempting to read with PetscBinaryRead...\n');
    b = PetscBinaryRead(inputFile);
    if ~isempty(b)
      success = true;
      fprintf('Read vector of length %d via PetscBinaryRead.\n', numel(b));
    end
  catch err
    warning('PetscBinaryRead failed: %s', err.message);
  end
end

fid = fopen(outputFile, 'w');
if fid == -1
  error('Could not open output file for writing: %s', outputFile);
end

n = numel(b);
for idx = 1:n
  fprintf(fid, '%d %.17g\n', idx, b(idx));   % 1-based index, high precision value
end
fclose(fid);

fprintf('Wrote %d lines to %s (index value, 1-based index).\n', n, outputFile);