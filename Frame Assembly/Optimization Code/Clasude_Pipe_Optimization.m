function tube_cutting_stock
% TUBE_CUTTING_STOCK  Exact min-cost cutting-stock solver for the frame tubes.
% Chooses how many 3.75 ft and 7.5 ft stock bars to buy, and which members
% get cut from which bar, to minimize total DOLLARS while cutting all members.
% Requires the Optimization Toolbox (intlinprog).
%
% Formulation: grouped-integer assignment form of the 1D cutting-stock problem
% (Gilmore & Gomory, Operations Research 9(6), 1961).
%   a(i,j) = integer >=0 : # pieces of distinct length i cut from bar j
%   u(j)   = binary      : 1 if bar j is purchased
%   min  sum_j cost(j)*u(j)
%   s.t. sum_j a(i,j) = d(i)                (demand: cut every member)
%        sum_i p_eff(i)*a(i,j) <= L(j)*u(j) (capacity: contents fit the bar)
%        u(j) <= u(j-1) within each type    (symmetry breaking)

% ================= INPUTS YOU OWN — TUNE THESE =========================
% --- Member table (distinct lengths, ft) and how many of each are needed
p = [5.000;4.960;3.917;3.061;3.054;2.846;2.780;2.667;2.560;2.529; ...
     2.441;2.206;2.133;2.122;1.944;1.939;1.640;1.608;1.562;1.476; ...
     1.422;1.312;1.164;1.050;0.669;0.652;0.578;0.304];        % ft
d = [2;1;2;2;2;2;2;2;2;2;2;2;1;2;2;1;1;1;4;1;1;1;2;6;2;2;2;2]; % count

% --- Stock catalog: nominal length (ft) and price ($/bar) for each type
Lnom  = [7.50 ; 3.75];      % nominal purchasable lengths, ft
price = [ 36 ;  16 ];     % <-- PUT SummitRacing per-bar prices here ($)

% --- Physical margins (your call, in feet)
kerf_ft  = 0.05/12;   % saw kerf removed per cut  (~0.05 in default)
face_ft  = 0.375/12;  % material faced off EACH bar end (~3/8 in default)

% --- How many candidate bars of each type to allow (upper bounds).
%     Total footage floor ~ 15 of the 7.5 ft bars; give generous slack.
Nbars = [18 ; 8];     % [n_7.5 ; n_3.75] candidate bars
% ======================================================================

% ---- Apply the derates -----------------------------------------------
p_eff = p + kerf_ft;              % each piece consumes one kerf (conservative)
L     = Lnom - 2*face_ft;         % usable length after facing both ends

% ---- Build the flat list of candidate bars ---------------------------
barL    = []; barCost = []; barType = [];
for t = 1:numel(Lnom)
    barL    = [barL   ; repmat(L(t),     Nbars(t),1)]; %#ok<AGROW>
    barCost = [barCost; repmat(price(t), Nbars(t),1)]; %#ok<AGROW>
    barType = [barType; repmat(t,        Nbars(t),1)]; %#ok<AGROW>
end
n = numel(p);       % distinct lengths
M = numel(barL);    % candidate bars
nVar = n*M + M;     % [ vec(a) ; u ]

idxA = @(i,j) (j-1)*n + i;   % column-major linear index of a(i,j)
idxU = @(j)   n*M + j;       % index of u(j)

% ---- Objective: cost only on the u(j) --------------------------------
f = zeros(nVar,1);
for j = 1:M, f(idxU(j)) = barCost(j); end

% ---- Demand equalities: sum_j a(i,j) = d(i) --------------------------
Aeq = zeros(n,nVar); beq = d;
for i = 1:n
    for j = 1:M, Aeq(i,idxA(i,j)) = 1; end
end

% ---- Capacity inequalities: sum_i p_eff(i)*a(i,j) - L(j)*u(j) <= 0 ----
Acap = zeros(M,nVar); bcap = zeros(M,1);
for j = 1:M
    for i = 1:n, Acap(j,idxA(i,j)) = p_eff(i); end
    Acap(j,idxU(j)) = -barL(j);
end

% ---- Symmetry breaking: u(j) <= u(j-1) within each type --------------
Asym = []; bsym = [];
for t = 1:numel(Lnom)
    js = find(barType==t);
    for k = 2:numel(js)
        row = zeros(1,nVar);
        row(idxU(js(k)))   =  1;
        row(idxU(js(k-1))) = -1;
        Asym = [Asym; row]; bsym = [bsym; 0]; %#ok<AGROW>
    end
end

A = [Acap; Asym]; b = [bcap; bsym];

% ---- Bounds ----------------------------------------------------------
lb = zeros(nVar,1);
ub = zeros(nVar,1);
for i = 1:n
    for j = 1:M
        ub(idxA(i,j)) = min(d(i), floor(barL(j)/p_eff(i))); % can't exceed need or fit
    end
end
for j = 1:M, ub(idxU(j)) = 1; end   % u binary

intcon  = 1:nVar;                    % everything integer
opts    = optimoptions('intlinprog','Display','iter');

[x,fval,exitflag] = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub,opts);

% ================= REPORT =============================================
if exitflag <= 0
    warning('No optimal solution found (exitflag=%d). Raise Nbars or check prices.',exitflag);
    return
end
a = reshape(x(1:n*M), n, M);
u = round(x(n*M+1:end));

fprintf('\n================ ORDER ================\n');
for t = 1:numel(Lnom)
    cnt = sum(u==1 & barType==t);
    fprintf('  %2d x %.2f ft bars\n', cnt, Lnom(t));
end
fprintf('  Total cost: $%.2f\n', fval);

purchased_ft = sum(u.*barL) + sum(u)*2*face_ft;   % add faced ends back for gross ft
gross_ft     = sum(u.*Lnom);
fprintf('  Gross purchased length: %.2f ft\n', gross_ft);
fprintf('  Member length required: %.3f ft\n', sum(p.*d));
fprintf('  Total waste (scrap+kerf+facing): %.3f ft\n', gross_ft - sum(p.*d));

fprintf('\n============ CUT LIST (per used bar) ============\n');
for j = 1:M
    if u(j)==1
        used = find(a(:,j)>0);
        pieces = [];
        for i = used(:)'
            pieces = [pieces; repmat(p(i),round(a(i,j)),1)]; %#ok<AGROW>
        end
        fprintf('  %.2f ft bar: [%s]  leftover %.3f ft\n', ...
            Lnom(barType(j)), strjoin(compose('%.3f',pieces),', '), ...
            barL(j) - sum(p_eff(used).*a(used,j)));
    end
end
end