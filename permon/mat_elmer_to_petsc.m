args = argv();
if numel(args) >= 1
  inputFile = args{1};
else
  inputFile = 'a.dat';
end
if numel(args) >= 2
  outputFile = args{2};
else
  outputFile = 'a.bin';
end

% open and read triplet-format file: "i j value" per line
fid = fopen(inputFile, 'r');
if fid == -1
  error('Could not open input file: %s', inputFile);
end
% allow comments starting with #, flexible whitespace
data = textscan(fid, '%f %f %f', 'CommentStyle', '#', 'MultipleDelimsAsOne', true);
fclose(fid);

if isempty(data) || isempty(data{1})
  error('No numeric data found in %s', inputFile);
end

i = data{1};
j = data{2};
v = data{3};

imax = max(i);
jmax = max(j);
n = max(imax, jmax);
A = sparse(double(i), double(j), v, n, n);

% write using PetscBinaryWrite
if exist('PetscBinaryWrite', 'file') || exist('PetscBinaryWrite', 'builtin')
  PetscBinaryWrite(outputFile, A);
  fprintf('Wrote PETSc binary matrix to %s\n', outputFile);
else
  error(['PetscBinaryWrite not found. Install or add the Octave PETSc interface ', ...
         'that provides PetscBinaryWrite to the path.']);
end
