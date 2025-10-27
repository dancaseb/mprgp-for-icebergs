args = argv();
if numel(args) >= 1
  inputFile = args{1};
else
  inputFile = 'linsys_b.dat';
end
if numel(args) >= 2
  outputFile = args{2};
else
  outputFile = 'linsys_b.bin';
end

% open and read double-format file: "i j value" per line
fid = fopen(inputFile, 'r');
if fid == -1
  error('Could not open input file: %s', inputFile);
end
% allow comments starting with #, flexible whitespace
data = textscan(fid, '%f %f', 'CommentStyle', '#', 'MultipleDelimsAsOne', true);
fclose(fid);

if isempty(data) || isempty(data{1})
  error('No numeric data found in %s', inputFile);
end

i = data{1};
v = data{2};

imax = max(i);

n = imax;
b = zeros(n,1);
b(i) = v;

% write using PetscBinaryWrite
if exist('PetscBinaryWrite', 'file') || exist('PetscBinaryWrite', 'builtin')
  PetscBinaryWrite(outputFile, b);
  fprintf('Wrote PETSc binary matrix to %s\n', outputFile);
else
  error(['PetscBinaryWrite not found. Install or add the Octave PETSc interface ', ...
         'that provides PetscBinaryWrite to the path.']);
end
