%% Build viridis256.csv from matplotlib's source
ref = 'v3.11.2';   
url = ['https://raw.githubusercontent.com/matplotlib/matplotlib/' ref '/lib/matplotlib/_cm_listed.py'];
txt = webread(url, weboptions('ContentType', 'text'));

i = strfind(txt, '_viridis_data = [');           % start of the viridis block
assert(isscalar(i), 'Expected exactly one _viridis_data block')
j = i - 1 + strfind(txt(i:end), ']]');           % end of the block
block = txt(i:j(1));

vals = str2double(regexp(block, '-?\d+(\.\d+)?([eE][-+]?\d+)?', 'match'));   % every number, in order
cmap = reshape(vals, 3, []).';                   % one row per colour: R G B

assert(isequal(size(cmap), [256 3]), 'Expected 256 x 3, got %d x %d', size(cmap))
assert(all(cmap(:) >= 0 & cmap(:) <= 1), 'Values should be in [0, 1]')

writematrix(cmap, 'viridis256.csv')
fprintf('Saved viridis256.csv from matplotlib %s\n', ref)