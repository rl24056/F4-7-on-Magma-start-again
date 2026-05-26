SetSeed(3);

print "============================================================";
print "F4(7) ADAPTIVE WALTON (2,3,13) SEARCH";
print "Use actual observed classes: involution traces 1,2 and order-3 traces -1,1";
print "Do not force a to have trace -2.";
print "No full Normalizer(G,<g>). No artificial GL correction.";
print "============================================================";

/* ============================================================
   PARAMETERS
   ============================================================ */

A_TRACE_TARGETS := [1, 2];       // observed genuine involution traces
B_TRACE_TARGETS := [-1, 1];      // observed genuine order-3 traces

REPS_PER_TRACE      := 10;
REP_ATTEMPTS        := 10000;

PRODUCT_TRIALS      := 150000;
PROGRESS_INTERVAL   := 5000;

ENUM_CAP            := 1092;
STOP_AFTER_FIRST    := true;

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
   If <a,b> grows above 1092, reject immediately.
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
   BUILD REPS BY TRACE
   ============================================================ */

function AllTraceBucketsFull(repsByTrace, wanted)
    for L in repsByTrace do
        if #L lt wanted then
            return false;
        end if;
    end for;

    return true;
end function;

function BuildPrimeOrderTraceReps(G, p, traceTargets, wantedPerTrace, maxAttempts)
    Id := Identity(G);

    repsByTrace := [ [ G | ] : i in [1..#traceTargets] ];
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

                for k in [1..#traceTargets] do
                    if tx eq traceTargets[k] and #repsByTrace[k] lt wantedPerTrace then
                        temp := repsByTrace[k];
                        AddIfNew(~temp, x);
                        repsByTrace[k] := temp;
                    end if;
                end for;
            end if;
        end if;

        if trial mod 1000 eq 0 then
            print "rep-build p =", p, "trial", trial;
            for k in [1..#traceTargets] do
                print "  trace", traceTargets[k], "count", #repsByTrace[k];
            end for;
        end if;

        if AllTraceBucketsFull(repsByTrace, wantedPerTrace) then
            print "All target reps for p =", p, "built by trial", trial;
            return repsByTrace, seenTraces;
        end if;
    end for;

    return repsByTrace, seenTraces;
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
        print "<a,b> exceeded 1092. Reject.";
        return false;
    end if;

    if #Hlist ne 1092 then
        print "<a,b> closed but wrong size. Reject.";
        return false;
    end if;

    print "SUCCESS: <a,b> has size 1092.";
    print "Now extract f inside H with g^f = g^10.";

    FCands := [ G | ];

    for x in Hlist do
        if IsExactOrder(x, 6, Id) then
            if x^-1 * g * x eq g^10 then
                Append(~FCands, x);
            end if;
        end if;
    end for;

    print "Number of f candidates with g^f=g^10:", #FCands;

    if #FCands eq 0 then
        print "No f found for this chosen g. Reject.";
        return false;
    end if;

    TargetFCands := [ x : x in FCands | Profile6(x) eq <1,-1,-2> ];

    if #TargetFCands gt 0 then
        f := TargetFCands[1];
        print "Using target V26-profile f.";
    else
        f := FCands[1];
        print "No target-profile f found; using first f.";
    end if;

    print "Chosen f profile [Tr(f),Tr(f^2),Tr(f^3)] =", Profile6(f);

    Blist := NormalFormList13by6(g, f);

    print "Normal-form count {g^i f^j} =", #Blist;

    if #Blist ne 78 then
        print "B is not clean 13:6. Reject.";
        return false;
    end if;

    print "Searching for external t in H with f^t=f^-1 and Bruhat witness.";

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
            print "Trace(a), Trace(b) =", TrS(a), TrS(b);
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

    print "H and B found, but no Walton t for this alignment.";
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
   STEP 2: BUILD REPS
   ============================================================ */

print "============================================================";
print "STEP 2: Build adaptive trace representatives";
print "Involutions: traces", A_TRACE_TARGETS;
print "Order-3: traces", B_TRACE_TARGETS;
print "============================================================";

ARepsByTrace, ASeen := BuildPrimeOrderTraceReps(
    G, 2, A_TRACE_TARGETS, REPS_PER_TRACE, REP_ATTEMPTS
);

BRepsByTrace, BSeen := BuildPrimeOrderTraceReps(
    G, 3, B_TRACE_TARGETS, REPS_PER_TRACE, REP_ATTEMPTS
);

print "------------------------------------------------------------";
print "Involution traces seen:";
for k in Keys(ASeen) do
    print k, ASeen[k];
end for;

print "Order-3 traces seen:";
for k in Keys(BSeen) do
    print k, BSeen[k];
end for;

APool := [ G | ];
ATraceLabel := [ Integers() | ];

for k in [1..#A_TRACE_TARGETS] do
    for x in ARepsByTrace[k] do
        Append(~APool, x);
        Append(~ATraceLabel, A_TRACE_TARGETS[k]);
    end for;
end for;

BPool := [ G | ];
BTraceLabel := [ Integers() | ];

for k in [1..#B_TRACE_TARGETS] do
    for x in BRepsByTrace[k] do
        Append(~BPool, x);
        Append(~BTraceLabel, B_TRACE_TARGETS[k]);
    end for;
end for;

print "#APool =", #APool;
print "#BPool =", #BPool;

if #APool eq 0 or #BPool eq 0 then
    print "Missing required reps. Script stops safely.";
else

    /* ========================================================
       STEP 3: RANDOM CONJUGATE SEARCH
       ======================================================== */

    print "============================================================";
    print "STEP 3: Random-conjugate (2,3,13) search";
    print "Trying actual observed classes, not forced trace -2.";
    print "============================================================";

    productHits := 0;
    finalHit := false;

    for trial in [1..PRODUCT_TRIALS] do
        ai := Random(1, #APool);
        bi := Random(1, #BPool);

        ca := Random(G);
        cb := Random(G);

        a := APool[ai]^ca;
        b := BPool[bi]^cb;

        g1 := a * b;

        if g1 ne IdG and g1^13 eq IdG and TrS(g1) eq 0 then
            productHits +:= 1;

            print "Candidate product hit", productHits, "at trial", trial;
            print "case g=a*b; Trace(a),Trace(b) =", TrS(a), TrS(b);

            finalHit := AnalysePairHit(G, a, b, g1, "g=a*b");

            if finalHit and STOP_AFTER_FIRST then
                break;
            end if;
        end if;

        g2 := b * a;

        if g2 ne IdG and g2^13 eq IdG and TrS(g2) eq 0 then
            productHits +:= 1;

            print "Candidate product hit", productHits, "at trial", trial;
            print "case g=b*a; Trace(a),Trace(b) =", TrS(a), TrS(b);

            finalHit := AnalysePairHit(G, a, b, g2, "g=b*a");

            if finalHit and STOP_AFTER_FIRST then
                break;
            end if;
        end if;

        if trial mod PROGRESS_INTERVAL eq 0 then
            print "trial", trial, "productHits", productHits;
        end if;
    end for;

    print "============================================================";
    print "SEARCH FINISHED";
    print "Product trials =", PRODUCT_TRIALS;
    print "Product hits tested =", productHits;
    print "Final Walton triple found?", finalHit;
    print "============================================================";

    if not finalHit then
        print "No final hit in this adaptive run.";
        print "This only means the random-conjugate search did not hit yet.";
        print "Next increase PRODUCT_TRIALS, or move this exact script to CERES/HPC.";
    end if;
end if;

print "DONE.";
