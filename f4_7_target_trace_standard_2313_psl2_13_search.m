
SetSeed(2);

print "============================================================";
print "F4(7) WALTON STANDARD (2,3,13) SEARCH -- TARGET TRACE CLASSES";
print "Fresh Magma Online script";
print "No full Normalizer. No full subgroup order unless subgroup closes at 1092.";
print "============================================================";

/* ============================================================
   PARAMETERS
   ============================================================ */

NUM_A_REPS        := 12;       // seed reps for involution trace -2
NUM_B_REPS        := 12;       // seed reps for order-3 trace -1
REP_ATTEMPTS      := 8000;     // attempts to build seed reps
PRODUCT_TRIALS    := 120000;   // main random-conjugate product tests
PROGRESS_INTERVAL := 5000;
ENUM_CAP          := 1092;     // PSL2(13) order
STOP_AFTER_FIRST  := true;

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

function InList(x, L)
    for y in L do
        if x eq y then
            return true;
        end if;
    end for;
    return false;
end function;

procedure AddIfNew(~L, x)
    if not InList(x, L) then
        Append(~L, x);
    end if;
end procedure;

function ActionExponentOnG(g, x)
    gx := x^-1 * g * x;
    for k in [1..12] do
        if gx eq g^k then
            return k;
        end if;
    end for;
    return 0;
end function;

function NormalFormList13by6(g, f)
    L := [ Parent(g) | ];

    for i in [0..12] do
        for j in [0..5] do
            x := g^i * f^j;
            if not InList(x, L) then
                Append(~L, x);
            end if;
        end for;
    end for;

    return L;
end function;

/* ============================================================
   SAFE ENUMERATION UP TO CAP
   This avoids trying to compute huge #<a,b>.
   If subgroup grows beyond 1092, reject.
   ============================================================ */

function EnumerateGeneratedSubgroupUpToCap(gens, cap)
    P := Parent(gens[1]);
    Id := Identity(P);

    usegens := [ P | ];
    for x in gens do
        Append(~usegens, x);
        Append(~usegens, x^-1);
    end for;

    L := [ P | Id ];
    idx := 1;

    while idx le #L do
        x := L[idx];

        for s in usegens do
            y := x * s;

            if not InList(y, L) then
                Append(~L, y);

                if #L gt cap then
                    return L, false;
                end if;
            end if;
        end for;

        idx +:= 1;
    end while;

    return L, true;
end function;

/* ============================================================
   BUILD TARGET TRACE REPS
   ============================================================ */

function BuildPrimeOrderTraceReps(G, p, trTarget, wanted, maxAttempts)
    Id := Identity(G);
    reps := [ G | ];

    seenTraces := AssociativeArray();

    for trial in [1..maxAttempts] do
        r := Random(G);
        o := Order(r);

        if o mod p eq 0 then
            x := r^(o div p);

            if x ne Id and x^p eq Id then
                tx := TrS(x);
                key := Sprint(tx);

                if IsDefined(seenTraces, key) then
                    seenTraces[key] +:= 1;
                else
                    seenTraces[key] := 1;
                end if;

                if tx eq trTarget then
                    AddIfNew(~reps, x);

                    if #reps ge wanted then
                        print "Built reps for order", p, "trace", trTarget, "at trial", trial;
                        return reps, seenTraces;
                    end if;
                end if;
            end if;
        end if;

        if trial mod 1000 eq 0 then
            print "rep-build p =", p, "trial", trial, "target reps =", #reps;
        end if;
    end for;

    return reps, seenTraces;
end function;

/* ============================================================
   ANALYSE A POSSIBLE (2,3,13) HIT
   ============================================================ */

function AnalysePairHit(G, a, b, g, label)
    Id := Identity(G);

    print "------------------------------------------------------------";
    print "PRODUCT HIT:", label;
    print "Order(a), Trace(a) =", Order(a), TrS(a);
    print "Order(b), Trace(b) =", Order(b), TrS(b);
    print "Order(g), Trace(g) =", Order(g), TrS(g);

    print "Enumerating <a,b> up to cap 1092 ...";
    Hlist, closed := EnumerateGeneratedSubgroupUpToCap([a,b], 1092);

    print "Enumeration closed?", closed;
    print "Enumerated size =", #Hlist;

    if not closed then
        print "<a,b> exceeded 1092, so this is not the target PSL2(13).";
        return false;
    end if;

    if #Hlist ne 1092 then
        print "<a,b> closed but has wrong size, not PSL2(13).";
        return false;
    end if;

    print "SUCCESS: <a,b> has size 1092, so this is a genuine PSL2(13) candidate.";

    print "Searching for f inside this 1092-element subgroup...";
    FCands := [ G | ];

    for x in Hlist do
        if IsExactOrder(x, 6, Id) then
            if x^-1 * g * x eq g^10 then
                Append(~FCands, x);
            end if;
        end if;
    end for;

    print "Number of f candidates with g^f = g^10:", #FCands;

    if #FCands eq 0 then
        print "No f found. Candidate subgroup has PSL2(13) size but wrong chosen g/alignment.";
        return false;
    end if;

    TargetFCands := [ x : x in FCands | Profile6(x) eq <1,-1,-2> ];

    if #TargetFCands gt 0 then
        f := TargetFCands[1];
        print "Using target-profile f.";
    else
        f := FCands[1];
        print "No target-profile f found; using first f.";
    end if;

    print "Chosen f profile [Tr(f),Tr(f^2),Tr(f^3)] =", Profile6(f);

    Blist := NormalFormList13by6(g, f);

    print "Normal-form count {g^i f^j} =", #Blist;

    if #Blist ne 78 then
        print "B=<g,f> is not clean 13:6. Rejecting this alignment.";
        return false;
    end if;

    print "Searching for Walton external involution t...";

    TCands := [ G | ];

    for t in Hlist do
        if IsExactOrder(t, 2, Id) then
            if t^-1 * f * t eq f^-1 then
                if not InList(t, Blist) then
                    Append(~TCands, t);
                end if;
            end if;
        end if;
    end for;

    print "External exact-inverting t candidates =", #TCands;

    for t in TCands do
        exp := ActionExponentOnG(g, t);

        witness := false;
        wi := -1;
        wj := -1;

        for i in [0..12] do
            for j in [0..5] do
                bb := g^i * f^j;
                q := bb * t;

                if q ne Id and q^3 eq Id then
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
            print "============================================================";
            print "FINAL WALTON TRIPLE FOUND INSIDE GENUINE F4(7)";
            print "============================================================";
            print "|H| = 1092";
            print "|B| normal forms =", #Blist;
            print "Order(g) =", Order(g);
            print "Order(f) =", Order(f);
            print "Order(t) =", Order(t);
            print "Trace(g) =", TrS(g);
            print "Trace(f), Trace(f^2), Trace(f^3) =", Profile6(f);
            print "Trace(t) =", TrS(t);
            print "Trace(g*t) =", TrS(g*t);
            print "f^-1*g*f = g^10 ?", f^-1 * g * f eq g^10;
            print "t^-1*f*t = f^-1 ?", t^-1 * f * t eq f^-1;
            print "t in B ?", InList(t, Blist);
            print "t normalises <g>? exponent 0 means no:", exp;
            print "Bruhat witness b=g^i f^j:";
            print "i,j =", wi, wj;
            print "Order((g^i f^j)*t) =", Order((g^wi * f^wj) * t);
            print "V26 trace profile [g,f,f^2,f^3,t,g*t] =",
                  [ TrS(g), TrS(f), TrS(f^2), TrS(f^3), TrS(t), TrS(g*t) ];
            print "============================================================";
            return true;
        end if;
    end for;

    print "PSL2(13) subgroup and B=13:6 found, but no Bruhat t found for this alignment.";
    return false;
end function;

/* ============================================================
   STEP 1: BUILD GENUINE F4(7)
   ============================================================ */

print "STEP 1: Build genuine F4(7)";
F := GF(7);
G := ChevalleyGroup("F", 4, F);
IdG := Identity(G);

print "Degree(G) =", Degree(G);
print "BaseRing(G) =", BaseRing(G);
print "Ngens(G) =", Ngens(G);

/* ============================================================
   STEP 2: BUILD TARGET TRACE CLASS REPRESENTATIVES
   ============================================================ */

print "============================================================";
print "STEP 2: Build target trace representatives";
print "Need order-2 trace -2 and order-3 trace -1.";
print "============================================================";

AReps, ATraceSeen := BuildPrimeOrderTraceReps(G, 2, -2, NUM_A_REPS, REP_ATTEMPTS);
BReps, BTraceSeen := BuildPrimeOrderTraceReps(G, 3, -1, NUM_B_REPS, REP_ATTEMPTS);

print "------------------------------------------------------------";
print "#AReps order-2 trace -2 =", #AReps;
print "#BReps order-3 trace -1 =", #BReps;

print "Involution traces seen while building AReps:";
for k in Keys(ATraceSeen) do
    print k, ATraceSeen[k];
end for;

print "Order-3 traces seen while building BReps:";
for k in Keys(BTraceSeen) do
    print k, BTraceSeen[k];
end for;

if #AReps eq 0 or #BReps eq 0 then
    print "Not enough target trace reps. Increase REP_ATTEMPTS or check trace convention.";
    error "Stopping: missing target trace class reps.";
end if;

/* ============================================================
   STEP 3: RANDOM CONJUGATE (2,3,13) PRODUCT SEARCH
   ============================================================ */

print "============================================================";
print "STEP 3: Random-conjugate standard-presentation search";
print "Search: |a|=2, Tr(a)=-2, |b|=3, Tr(b)=-1, |ab|=13, Tr(ab)=0.";
print "============================================================";

productHits := 0;
finalHit := false;

for trial in [1..PRODUCT_TRIALS] do
    ai := Random(1, #AReps);
    bi := Random(1, #BReps);

    ca := Random(G);
    cb := Random(G);

    a := AReps[ai]^ca;
    b := BReps[bi]^cb;

    g1 := a * b;

    if g1 ne IdG and g1^13 eq IdG and TrS(g1) eq 0 then
        productHits +:= 1;
        print "Candidate product hit number", productHits, "at trial", trial, "for g=a*b";

        finalHit := AnalysePairHit(G, a, b, g1, "g=a*b");

        if finalHit and STOP_AFTER_FIRST then
            break;
        end if;
    end if;

    g2 := b * a;

    if g2 ne IdG and g2^13 eq IdG and TrS(g2) eq 0 then
        productHits +:= 1;
        print "Candidate product hit number", productHits, "at trial", trial, "for g=b*a";

        finalHit := AnalysePairHit(G, a, b, g2, "g=b*a");

        if finalHit and STOP_AFTER_FIRST then
            break;
        end if;
    end if;

    if trial mod PROGRESS_INTERVAL eq 0 then
        print "trial", trial, "productHits so far", productHits;
    end if;
end for;

print "================================================------------";
print "SEARCH FINISHED";
print "Product trials =", PRODUCT_TRIALS;
print "Product hits tested =", productHits;
print "Final Walton triple found?", finalHit;
print "================================================------------";

if not finalHit then
    print "No final hit in this run.";
    print "Next controlled increases:";
    print "1. PRODUCT_TRIALS := 300000";
    print "2. NUM_A_REPS := 20; NUM_B_REPS := 20";
    print "3. REP_ATTEMPTS := 20000";
    print "Do not use Normalizer(G,<g>) in Magma Online.";
end if;

print "DONE.";
