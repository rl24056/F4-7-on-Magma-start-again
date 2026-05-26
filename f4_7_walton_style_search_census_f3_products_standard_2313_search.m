SetSeed(1);

print "============================================================";
print "F4(7) WALTON-STYLE SEARCH: CENSUS + f^3 PRODUCTS + (2,3,13)";
print "No full Normalizer in F4(7). No artificial GL/Hermitian correction.";
print "============================================================";

/* ============================================================
   PARAMETERS
   ============================================================ */

POOL_TRIALS      := 3000;   // random elements used to build pools
MAX_INV_POOL     := 80;     // involution pool cap
MAX_ORD3_POOL    := 80;     // order-3 pool cap
MAX_ORD6_POOL    := 80;     // order-6 pool cap
MAX_PAIR_TESTS   := 8000;   // cap for (2,3,13) product tests
MAX_F3_TESTS     := 8000;   // cap for f^3 involution product tests

/* ============================================================
   BASIC HELPERS
   ============================================================ */

function SignedGF7(x)
    z := Integers()!x;
    if z ge 4 then
        return z - 7;
    else
        return z;
    end if;
end function;

function TrS(x)
    return SignedGF7(Trace(Matrix(x)));
end function;

function Profile6(x)
    return <TrS(x), TrS(x^2), TrS(x^3)>;
end function;

function IsExactOrder(x, n, Id)
    if x eq Id then
        return false;
    end if;
    if x^n ne Id then
        return false;
    end if;
    for p in PrimeDivisors(n) do
        if x^(n div p) eq Id then
            return false;
        end if;
    end for;
    return true;
end function;

procedure AddIfNew(~S, x, cap)
    if #S ge cap then
        return;
    end if;

    for y in S do
        if x eq y then
            return;
        end if;
    end for;

    Append(~S, x);
end procedure;

function NormalFormCount13by6(g, f)
    S := [ Parent(g) | ];
    for i in [0..12] do
        for j in [0..5] do
            x := g^i * f^j;

            seen := false;
            for y in S do
                if x eq y then
                    seen := true;
                    break;
                end if;
            end for;

            if not seen then
                Append(~S, x);
            end if;
        end for;
    end for;

    return #S;
end function;

function ActionExponentOnG(g, x, Id)
    gx := x^-1 * g * x;
    for k in [1..12] do
        if gx eq g^k then
            return k;
        end if;
    end for;
    return 0;
end function;

/* ============================================================
   ANALYSE A SMALL STANDARD-PRESENTATION HIT
   a^2 = b^3 = (ab)^13 = 1
   ============================================================ */

function AnalyseStandardHit(G, a, b, g, label)
    print "------------------------------------------------------------";
    print "STANDARD HIT FOUND:", label;
    print "Trace(g) =", TrS(g);
    print "Checking H = <a,b> ...";

    H := sub< G | a, b >;
    Hsize := #H;

    print "#H =", Hsize;

    if Hsize ne 1092 then
        print "H is not PSL2(13) of order 1092. Interesting, but not final target.";
        return false;
    end if;

    print "SUCCESS: #H = 1092, so this is the correct PSL2(13) size.";

    gH := H!g;
    IdH := Identity(H);
    EH := Elements(H);

    print "Searching for f inside H with order 6 and g^f = g^10 ...";

    FCands := [ H | ];

    for x in EH do
        if IsExactOrder(x, 6, IdH) then
            if x^-1 * gH * x eq gH^10 then
                Append(~FCands, x);
            end if;
        end if;
    end for;

    print "Number of f candidates in H =", #FCands;

    if #FCands eq 0 then
        print "No f found inside H. This should not happen if the PSL2(13) alignment is right.";
        return false;
    end if;

    TargetFCands := [ x : x in FCands | Profile6(x) eq <1,-1,-2> ];

    if #TargetFCands gt 0 then
        fH := TargetFCands[1];
        print "Using target-profile f with profile [Tr(f),Tr(f^2),Tr(f^3)] =", Profile6(fH);
    else
        fH := FCands[1];
        print "No target-profile f found; using first f.";
        print "Chosen f profile [Tr(f),Tr(f^2),Tr(f^3)] =", Profile6(fH);
    end if;

    B := sub< H | gH, fH >;
    print "#B = #<g,f> =", #B;
    print "Normal-form count {g^i f^j} =", NormalFormCount13by6(gH, fH);

    print "Searching for external t in H with f^t=f^-1 and Bruhat witness ...";

    TCands := [ H | ];

    for t in EH do
        if IsExactOrder(t, 2, IdH) then
            if t^-1 * fH * t eq fH^-1 then
                if not (t in B) then
                    Append(~TCands, t);
                end if;
            end if;
        end if;
    end for;

    print "External exact-inverting t candidates =", #TCands;

    for tH in TCands do
        exp := ActionExponentOnG(gH, tH, IdH);

        witness := false;
        wi := -1;
        wj := -1;

        for i in [0..12] do
            for j in [0..5] do
                bb := gH^i * fH^j;
                q := bb * tH;

                if q ne IdH and q^3 eq IdH then
                    witness := true;
                    wi := i;
                    wj := j;
                    break;
                end if;
            end for;

            if witness then
                break;
            end if;
        end for;

        if witness then
            print "------------------------------------------------------------";
            print "FINAL WALTON TRIPLE FOUND INSIDE GENUINE F4(7)";
            print "#H =", #H;
            print "#B =", #B;
            print "Order(g) =", Order(gH);
            print "Order(f) =", Order(fH);
            print "Order(t) =", Order(tH);
            print "f^-1*g*f = g^10 ?", fH^-1 * gH * fH eq gH^10;
            print "f^t = f^-1 ?", tH^-1 * fH * tH eq fH^-1;
            print "t in B ?", tH in B;
            print "t normalises <g>? exponent 0 means no:", exp;
            print "Bruhat witness: b = g^i f^j with (b*t)^3=1";
            print "i,j =", wi, wj;
            print "Order(g*t) =", Order(gH * tH);
            print "V26 trace profile [g,f,f^2,f^3,t,g*t] =",
                  [ TrS(gH), TrS(fH), TrS(fH^2), TrS(fH^3), TrS(tH), TrS(gH*tH) ];
            print "------------------------------------------------------------";
            return true;
        end if;
    end for;

    print "H=PSL2(13) found, f found, but no Bruhat t found in this alignment.";
    return false;
end function;

/* ============================================================
   ANALYSE LOCAL 13:6 HIT FROM f^3 PRODUCT ROUTE
   ============================================================ */

function AnalyseLocal13by6(G, g, f, label)
    IdG := Identity(G);

    if not IsExactOrder(g, 13, IdG) then
        return false;
    end if;

    if not IsExactOrder(f, 6, IdG) then
        return false;
    end if;

    if f^-1 * g * f ne g^10 then
        return false;
    end if;

    nf := NormalFormCount13by6(g, f);

    print "------------------------------------------------------------";
    print "LOCAL 13:6 HIT:", label;
    print "Order(g) =", Order(g);
    print "Order(f) =", Order(f);
    print "Trace(g) =", TrS(g);
    print "f^-1*g*f = g^10 ?", f^-1 * g * f eq g^10;
    print "f profile [Tr(f),Tr(f^2),Tr(f^3)] =", Profile6(f);
    print "Normal-form count {g^i f^j} =", nf;
    print "------------------------------------------------------------";

    if nf eq 78 then
        print "SUCCESS: genuine local B=<g,f> has 78 normal forms.";
        return true;
    else
        print "Not a clean 13:6 normal-form subgroup.";
        return false;
    end if;
end function;

/* ============================================================
   STEP 1: BUILD GENUINE F4(7)
   ============================================================ */

print "STEP 1: Building genuine F4(7) on V26...";
F := GF(7);
G := ChevalleyGroup("F", 4, F);
IdG := Identity(G);

print "Degree(G) =", Degree(G);
print "Base ring =", BaseRing(G);
print "Ngens(G) =", Ngens(G);

/* ============================================================
   STEP 2: BUILD POOLS OF GENUINE ORDER 2, 3, 6 ELEMENTS
   ============================================================ */

print "============================================================";
print "STEP 2: Random genuine class census / pool construction";
print "============================================================";

InvPool := [ G | ];
Ord3Pool := [ G | ];
Ord6Pool := [ G | ];

ProfileCount := AssociativeArray();

for trial in [1..POOL_TRIALS] do
    r := Random(G);
    o := Order(r);

    if o mod 2 eq 0 then
        a := r^(o div 2);
        if IsExactOrder(a, 2, IdG) then
            AddIfNew(~InvPool, a, MAX_INV_POOL);
        end if;
    end if;

    if o mod 3 eq 0 then
        b := r^(o div 3);
        if IsExactOrder(b, 3, IdG) then
            AddIfNew(~Ord3Pool, b, MAX_ORD3_POOL);
        end if;
    end if;

    if o mod 6 eq 0 then
        f := r^(o div 6);
        if IsExactOrder(f, 6, IdG) then
            AddIfNew(~Ord6Pool, f, MAX_ORD6_POOL);

            key := Sprint(Profile6(f));
            if IsDefined(ProfileCount, key) then
                ProfileCount[key] +:= 1;
            else
                ProfileCount[key] := 1;
            end if;
        end if;
    end if;

    if trial mod 250 eq 0 then
        printf "trial %o: involutions=%o, order3=%o, order6=%o\n",
               trial, #InvPool, #Ord3Pool, #Ord6Pool;
    end if;

    if #InvPool ge MAX_INV_POOL and #Ord3Pool ge MAX_ORD3_POOL and #Ord6Pool ge MAX_ORD6_POOL then
        print "All pools reached cap early.";
        break;
    end if;
end for;

print "------------------------------------------------------------";
print "Pool summary:";
print "#InvPool =", #InvPool;
print "#Ord3Pool =", #Ord3Pool;
print "#Ord6Pool =", #Ord6Pool;

print "Order-6 trace-profile counts from random census:";
for k in Keys(ProfileCount) do
    print k, ProfileCount[k];
end for;

/* ============================================================
   STEP 3: f^3 INVOLUTION PRODUCT ROUTE
   y=f^3, z another involution, g=z*y or y*z
   Test order(g)=13, Trace(g)=0, and f^-1*g*f=g^10
   ============================================================ */

print "============================================================";
print "STEP 3: f^3 involution-product route";
print "============================================================";

f3Tests := 0;
localHit := false;

for i in [1..#Ord6Pool] do
    f := Ord6Pool[i];
    y := f^3;

    for j in [1..#InvPool] do
        z := InvPool[j];

        if z eq y then
            continue;
        end if;

        for side in [1..2] do
            if side eq 1 then
                g := z * y;
                label := "g = z*f^3";
            else
                g := y * z;
                label := "g = f^3*z";
            end if;

            f3Tests +:= 1;

            if IsExactOrder(g, 13, IdG) and TrS(g) eq 0 then
                print "Candidate from f^3 route:", label;
                print "f profile =", Profile6(f);
                print "Trace(g) =", TrS(g);

                if f^-1 * g * f eq g^10 then
                    localHit := AnalyseLocal13by6(G, g, f, label);
                    if localHit then
                        break;
                    end if;
                else
                    print "But f does not act by exponent 10 on g.";
                end if;
            end if;

            if f3Tests ge MAX_F3_TESTS then
                break;
            end if;
        end for;

        if localHit or f3Tests ge MAX_F3_TESTS then
            break;
        end if;
    end for;

    if localHit or f3Tests ge MAX_F3_TESTS then
        break;
    end if;
end for;

print "f^3 route tests completed =", f3Tests;

if not localHit then
    print "No clean local 13:6 found from f^3 route in this capped run.";
end if;

/* ============================================================
   STEP 4: WALTON STANDARD-PRESENTATION ROUTE
   Search a,b with |a|=2, |b|=3, |ab|=13, Trace(ab)=0.
   If hit, analyse H=<a,b> and extract f,t inside H.
   ============================================================ */

print "============================================================";
print "STEP 4: Standard (2,3,13) product search";
print "============================================================";

if #InvPool eq 0 or #Ord3Pool eq 0 then
    print "Insufficient pools for (2,3,13) search.";
else
    pairTests := 0;
    finalHit := false;

    for ia in [1..#InvPool] do
        a := InvPool[ia];

        for ib in [1..#Ord3Pool] do
            b := Ord3Pool[ib];

            pairTests +:= 1;

            g1 := a * b;

            if IsExactOrder(g1, 13, IdG) and TrS(g1) eq 0 then
                print "Candidate standard pair: g = a*b";
                print "ia, ib =", ia, ib;
                finalHit := AnalyseStandardHit(G, a, b, g1, "g=a*b");

                if finalHit then
                    break;
                end if;
            end if;

            g2 := b * a;

            if IsExactOrder(g2, 13, IdG) and TrS(g2) eq 0 then
                print "Candidate standard pair: g = b*a";
                print "ia, ib =", ia, ib;
                finalHit := AnalyseStandardHit(G, a, b, g2, "g=b*a");

                if finalHit then
                    break;
                end if;
            end if;

            if pairTests ge MAX_PAIR_TESTS then
                break;
            end if;
        end for;

        if finalHit or pairTests ge MAX_PAIR_TESTS then
            break;
        end if;
    end for;

    print "Standard-pair tests completed =", pairTests;

    if finalHit then
        print "============================================================";
        print "FINAL SUCCESS IN THIS RUN.";
        print "A genuine PSL2(13) Walton triple was found inside F4(7).";
        print "============================================================";
    else
        print "============================================================";
        print "No final PSL2(13) Walton triple found in this capped run.";
        print "Next increase POOL_TRIALS, MAX_INV_POOL, MAX_ORD3_POOL, MAX_PAIR_TESTS.";
        print "Do NOT switch back to full Normalizer(G,<g>) in Magma Online.";
        print "============================================================";
    end if;
end if;

print "DONE.";
