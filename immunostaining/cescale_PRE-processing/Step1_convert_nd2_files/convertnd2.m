%% Step 1: ND2 -> hyb0_channels/*.mat
%Writes convert_nd2_v2.py and runs it under a module-loaded python.
%Edit THIS file, not the .py -- the .py is regenerated on every run.

%% ---- EDIT THESE  ----
MODULES = "python/3.11 scipy-stack";   %% python/3.11 is the one that has nd2
ND2     = "/path/to/nd2"; %your full path to the .nd2 file
WORKDIR = "/path/to/workdir";   %the directory of your .nd2 file
DRY_RUN = "False";    %% "True" = print metadata and write nothing. "False" = convert
FOVS    = [6];   %% which FOVs to convert; [] = all.
%% --------------------------

RUN = @(c) system(char(sprintf('bash -lc ''module load %s && %s''', MODULES, c)));

%% ---- pre-run ----
[st, out] = RUN('python -c "import nd2, scipy; print(nd2.__version__, scipy.__version__)"');
if st ~= 0
    fprintf('%s\n', out);
    error(['Preflight failed under modules "%s". If scipy is the problem, drop ' ...
           'scipy-stack from MODULES and install it offline (no login node needed):\n' ...
           '   bash -lc ''module load python/3.11 && pip install --no-index --user scipy'''], MODULES);
end
fprintf('nd2 / scipy versions: %s\n', strtrim(out));

%% ---- write the python file ----
pylines = [ ...
"import nd2"
"import numpy as np"
"import os"
"import re"
"from scipy.io import savemat"
""
"ND2     = '__ND2__'"
"WORKDIR = '__OUT__'"
"DRY_RUN = __DRY__"
"FOVS    = __FOVS__      # [] = every FOV in the file"
"HYB     = 0"
""
"ROLES = {'405': ('DAPI',     'ImageStackDAPI'),"
"         '488': ('Lineage',  'ImageStackLineage'),"
"         '561': ('PGranule', 'ImageStackPGranule'),"
"         '640': ('PGranule', 'ImageStackPGranule'),"
"         '647': ('PGranule', 'ImageStackPGranule')}"
""
"f = nd2.ND2File(ND2)"
"print('file   :', os.path.basename(ND2))"
"print('sizes  :', dict(f.sizes))"
"print('dtype  :', f.dtype)"
"def chname(c):"
"    try:"
"        return c.channel.name"
"    except AttributeError:"
"        return str(c)"
"names = [chname(c) for c in f.metadata.channels]"
"print('channels (file order):', names)"
"print('raw channel meta   :', f.metadata.channels)"
"try:"
"    print('voxel size (um, x y z):', f.voxel_size())"
"except Exception as e:"
"    print('voxel size unavailable:', e)"
""
"dims = list(f.sizes.keys())"
"nP   = f.sizes.get('P', 1)"
"nC   = f.sizes['C']"
"assert 'C' in dims and 'Z' in dims, dims"
"assert len(names) == nC, (names, nC)"
""
"roles = []"
"for n in names:"
"    hit = [w for w in ROLES if w in re.sub(r'[^0-9]', ' ', n)]"
"    if len(hit) != 1:"
"        raise SystemExit('cannot map channel %r to a role (matched %r)' % (n, hit))"
"    roles.append(ROLES[hit[0]])"
"for n, (p, v) in zip(names, roles):"
"    print('  %-12s -> %s_%d_{N}.mat   %s' % (n, p, HYB, v))"
"assert len({p for p, _ in roles}) == nC, 'two channels mapped to the same role'"
"want = list(range(1, nP + 1)) if not FOVS else sorted(set(FOVS))"
"bad  = [n for n in want if n < 1 or n > nP]"
"if bad:"
"    raise SystemExit('FOVS out of range 1..%d: %s' % (nP, bad))"
"print('FOVs   : %d in file, converting %d: %s' % (nP, len(want), want))"
""
"if DRY_RUN:"
"    print('DRY_RUN -- nothing written.')"
"    f.close()"
"    raise SystemExit"
""
"outdir = os.path.join(WORKDIR, 'hyb0_channels')"
"os.makedirs(outdir, exist_ok=True)"
"d = f.to_dask()"
""
"for n in want:"
"    p = n - 1"
"    for c in range(nC):"
"        prefix, varname = roles[c]"
"        path = os.path.join(outdir, '%s_%d_%d.mat' % (prefix, HYB, n))"
"        if os.path.exists(path):"
"            print('skip (exists):', os.path.basename(path), flush=True)"
"            continue"
"        idx = tuple(p if dim == 'P' else c if dim == 'C' else slice(None) for dim in dims)"
"        stack = np.asarray(d[idx].compute())"
"        assert stack.ndim == 3, (stack.shape, dims)"
"        stack = np.ascontiguousarray(np.transpose(stack, (1, 2, 0)))"
"        if stack.dtype != np.uint16:"
"            print('WARNING: %s is %s, casting to uint16' % (names[c], stack.dtype))"
"            stack = stack.astype(np.uint16)"
"        if stack.nbytes > 3.5e9:"
"            raise SystemExit('stack is %.1f GB, too big for MAT v5; needs an HDF5 writer'"
"                             % (stack.nbytes / 1e9))"
"        savemat(path, {varname: stack}, do_compression=False)"
"        print('%-9s N=%-3d  %s  %s  %.2f GB' % (prefix, n, stack.shape,"
"                                                stack.dtype, stack.nbytes / 1e9), flush=True)"
""
"f.close()"
"print('Done')"
];

script = strjoin(pylines, newline);
script = strrep(script, "__ND2__", ND2);
script = strrep(script, "__OUT__", WORKDIR);
script = strrep(script, "__DRY__", DRY_RUN);
if isempty(FOVS)
    fovstr = "[]";
else
    fovstr = "[" + strjoin(string(FOVS(:)'), ", ") + "]";
end
script = strrep(script, "__FOVS__", fovstr);

if contains(script, "__ND2__") || contains(script, "__OUT__") || contains(script, "__DRY__") || contains(script, "__FOVS__")
    error('Placeholder substitution failed -- the embedded python still has __ND2__/__OUT__/__DRY__/__FOVS__.');
end

fid = fopen('convert_nd2_v2.py', 'w');
fprintf(fid, '%s', script);
fclose(fid);

%% ---- run it ----
status = RUN('python convert_nd2_v2.py');
if status ~= 0
    error('convert_nd2_v2.py exited with status %d', status);
end