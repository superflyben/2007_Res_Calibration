function Build_2000_Surface(method)
%Build the 2000 Juan de Fuca calibration surface and save it for the R scripts.
%
%   Build_2000_Surface          % cubic (the default; see below)
%   Build_2000_Surface('linear')
%
%Derived from TSfit3.m (Udub/Modeling/3rdTSfit/), which produced the surface
%behind the GCA paper. The original machinery is kept deliberately: the five
%isohaline files are read and sorted by temperature as readlarson.m did, and
%interpmmlarson.m splines conductivity against temperature within each isohaline
%out to 10000 points. That pre-pass is what carries the shape of C(T), which
%rises to a plateau near 300 C and falls away sharply above ~350 C -- four
%measured temperatures per isohaline cannot represent that any other way.
%
%ONE deliberate change from TSfit3.m: gridfit is replaced by griddata. gridfit
%extrapolates, leaving the archived surface defined across its whole rectangle
%with no NaN, including regions with no data behind them. griddata instead
%returns NaN outside the convex hull of the data, which is what lets the
%comparison against the 2007 surface be masked to real measurements.
%
%Method: linear and cubic are near-indistinguishable here -- median relative
%difference 0.05%, max 1.3%, identical fidelity to the 20 measured points, and
%neither overshoots the data range. The spline pre-pass leaves the point cloud
%dense enough along each isohaline that the gridding method hardly matters.
%Cubic is the default, for the smooth surface the original fit was after.
%
%Salt is converted from mol/kg to mmol/kg on ingest so the output matches the
%2007 data and the R scripts downstream.
%
%Run with the working directory set to this folder; it reads bare filenames.

if nargin < 1
    method = 'cubic';
end

%--- read the constant-salt lines (readlarson.m, inlined so paths stay local)
datadir = fullfile('..','data','2000_JdF');
names = {'mm54','mm111','mm272','mm463','mm634'};
alldat = [];
numinterp = 10000;
for k = 1:numel(names)
    iso = csvread(fullfile(datadir, [names{k} '.csv']));
    iso = sortrows(iso, 2);
    %spline C against T within this isohaline
    alldat = [alldat; interpmmlarson(iso, numinterp)];
end

S = alldat(:,1) * 1000;   %concentration, mol/kg -> mmol/kg
T = alldat(:,2);          %temperature, deg C
C = alldat(:,3);          %conductivity, 1/V

%--- the grid of record: the nodes the 2007 Matlab surface uses, so the two
%    calibrations sit on identical nodes and can be differenced directly.
%    Verified to reproduce Calibration/2007/surfdat.mat's grid exactly.
[Tdat, Sdat] = meshgrid(linspace(202, 395.5, 20), linspace(100, 1500, 20));
Cdat = griddata(T, S, C, Tdat, Sdat, method);

%--- a dense surface over the data's own range, for plotting
[Tnat, Snat] = meshgrid(linspace(min(T), max(T), 200), ...
                        linspace(min(S), max(S), 200));
Cnat = griddata(T, S, C, Tnat, Snat, method);

%--- the measured points, carried along so R can plot them on the surface
Traw = zeros(0,1); Sraw = zeros(0,1); Craw = zeros(0,1);
for k = 1:numel(names)
    iso = sortrows(csvread(fullfile(datadir, [names{k} '.csv'])), 2);
    Sraw = [Sraw; iso(:,1)*1000]; Traw = [Traw; iso(:,2)]; Craw = [Craw; iso(:,3)];
end

save('surfdat2000.mat', 'Tdat','Sdat','Cdat','Tnat','Snat','Cnat', ...
     'Traw','Sraw','Craw','method');

fprintf('Built the 2000 surface with griddata ''%s''\n', method);
fprintf('  grid of record: %d of %d nodes finite\n', ...
    sum(isfinite(Cdat(:))), numel(Cdat));
fprintf('  native 200x200: %d of %d nodes finite\n', ...
    sum(isfinite(Cnat(:))), numel(Cnat));
fprintf('  C over the data: %.4f to %.4f; over the surface: %.4f to %.4f\n', ...
    min(Craw), max(Craw), min(Cdat(:)), max(Cdat(:)));
fprintf('  wrote surfdat2000.mat\n');
end
