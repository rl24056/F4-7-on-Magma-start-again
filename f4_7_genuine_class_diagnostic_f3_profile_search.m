SetSeed(7);

print "============================================================";
print "F4(7) GENUINE CLASS DIAGNOSTIC + f^3 PRODUCT SEARCH BY PROFILE";
print "No Normalizer. No PSL2 abstract code. Magma Online safe.";
print "============================================================";

/* parameters */
RANDOM_TRIALS := 10000;
MAX_REPS_PER_PROFILE := 3;
F3_TRIALS_PER_F := 2000;
PROGRESS := 1000;

/* helpers */

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

function FixDim(x)
    M := Matrix(x);
    F := BaseRing(M);
    n := Nrows(M);
    I := IdentityMatrix(F,n);
    return n - Rank(M - I);
end function;

function IsExactOrder(x,n,Id)
    if x eq Id then return false; end if;
    if x^n ne Id then return false; end if;

    for p in PrimeDivisors(n) do
        if x^(n div p) eq Id then
            return false;
        end if;
    end for;

    return true;
end function;

function InList(x,L)
    for y in L do
        if x eq y then return true; end if;
    end for;
    return false;
end function;

procedure AddIfNew(~L,x)
    if not InList(x,L) then
        Append(~L,x);
    end if;
end procedure;

function Prof6(f)
    return <TrS(f), TrS(f^2), TrS(f^3), FixDim(f), FixDim(f^2), FixDim(f^3)>;
end function;

function ActionExponent(g,x)
    gx := x^-1*g*x;

    for k in [1..12] do
        if gx eq g^k then
            return k;
        end if;
    end for;

    return 0;
end function;

function NormalFormCount(g,f)
    L := [ Parent(g) | ];

    for i in [0..12] do
        for j in [0..5] do
            x := g^i*f^j;

            if not InList(x,L) then
                Append(~L,x);
            end if;
        end for;
    end for;

    return #L;
end function;

/* build genuine F4(7) */

F := GF(7);
G := ChevalleyGroup("F",4,F);
Id := Identity(G);

print "Degree(G) =", Degree(G);
print "BaseRing(G) =", BaseRing(G);
print "Ngens(G) =", Ngens(G);

/* census */

print "============================================================";
print "STEP 1: random class/profile census";
print "============================================================";

InvTraceCount := AssociativeArray();
Ord3TraceCount := AssociativeArray();
Ord13TraceCount := AssociativeArray();
Ord6ProfileCount := AssociativeArray();
Ord6Reps := AssociativeArray();

for trial in [1..RANDOM_TRIALS] do
    r := Random(G);
    o := Order(r);

    if o mod 2 eq 0 then
        a := r^(o div 2);
        if IsExactOrder(a,2,Id) then
            k := Sprint(<TrS(a),FixDim(a)>);
            if IsDefined(InvTraceCount,k) then
                InvTraceCount[k] +:= 1;
            else
                InvTraceCount[k] := 1;
            end if;
        end if;
    end if;

    if o mod 3 eq 0 then
        b := r^(o div 3);
        if IsExactOrder(b,3,Id) then
            k := Sprint(<TrS(b),FixDim(b)>);
            if IsDefined(Ord3TraceCount,k) then
                Ord3TraceCount[k] +:= 1;
            else
                Ord3TraceCount[k] := 1;
            end if;
        end if;
    end if;

    if o mod 13 eq 0 then
        g := r^(o div 13);
        if IsExactOrder(g,13,Id) then
            k := Sprint(<TrS(g),FixDim(g)>);
            if IsDefined(Ord13TraceCount,k) then
                Ord13TraceCount[k] +:= 1;
            else
                Ord13TraceCount[k] := 1;
            end if;
        end if;
    end if;

    if o mod 6 eq 0 then
        f := r^(o div 6);

        if IsExactOrder(f,6,Id) then
            p := Prof6(f);
            k := Sprint(p);

            if IsDefined(Ord6ProfileCount,k) then
                Ord6ProfileCount[k] +:= 1;
            else
                Ord6ProfileCount[k] := 1;
            end if;

            if IsDefined(Ord6Reps,k) then
                L := Ord6Reps[k];
                if #L lt MAX_REPS_PER_PROFILE then
                    AddIfNew(~L,f);
                    Ord6Reps[k] := L;
                end if;
            else
                Ord6Reps[k] := [ G | f ];
            end if;
        end if;
    end if;

    if trial mod PROGRESS eq 0 then
        print "trial", trial;
    end if;
end for;

print "------------------------------------------------------------";
print "Involution profiles <trace,fixdim>:";
for k in Keys(InvTraceCount) do
    print k, InvTraceCount[k];
end for;

print "------------------------------------------------------------";
print "Order-3 profiles <trace,fixdim>:";
for k in Keys(Ord3TraceCount) do
    print k, Ord3TraceCount[k];
end for;

print "------------------------------------------------------------";
print "Order-13 profiles <trace,fixdim>:";
for k in Keys(Ord13TraceCount) do
    print k, Ord13TraceCount[k];
end for;

print "------------------------------------------------------------";
print "Order-6 profiles <Tr(f),Tr(f^2),Tr(f^3),Fix(f),Fix(f^2),Fix(f^3)>:";
for k in Keys(Ord6ProfileCount) do
    print k, Ord6ProfileCount[k];
end for;

print "============================================================";
print "STEP 2: f^3 product search by genuine order-6 profile";
print "For each f, set y=f^3, z=y^c, g=z*y or y*z.";
print "Test order(g)=13, Tr(g)=0, and whether f normalises <g>.";
print "============================================================";

totalProductHits := 0;
normalisingHits := 0;
localBHits := 0;

for key in Keys(Ord6Reps) do
    print "------------------------------------------------------------";
    print "Testing order-6 profile:", key;

    reps := Ord6Reps[key];

    for rr in [1..#reps] do
        f := reps[rr];
        y := f^3;

        print "  rep", rr, "Trace(y=f^3), Fix(y) =", TrS(y), FixDim(y);

        profileProductHits := 0;
        profileNormalisingHits := 0;

        for trial in [1..F3_TRIALS_PER_F] do
            c := Random(G);
            z := y^c;

            for side in [1..2] do
                if side eq 1 then
                    g := z*y;
                    label := "g=z*y";
                else
                    g := y*z;
                    label := "g=y*z";
                end if;

                if g ne Id and g^13 eq Id and TrS(g) eq 0 then
                    totalProductHits +:= 1;
                    profileProductHits +:= 1;

                    exp := ActionExponent(g,f);

                    if exp ne 0 then
                        normalisingHits +:= 1;
                        profileNormalisingHits +:= 1;

                        print "  PRODUCT+NORMALISING HIT";
                        print "  label =", label;
                        print "  trial =", trial;
                        print "  exp of f on <g> =", exp;
                        print "  Trace(g), Fix(g) =", TrS(g), FixDim(g);
                        print "  f profile =", Prof6(f);

                        if exp eq 10 then
                            ff := f;
                            nf := NormalFormCount(g,ff);

                            print "  Normal-form count for <g,f> =", nf;

                            if nf eq 78 then
                                localBHits +:= 1;
                                print "============================================================";
                                print "LOCAL B=13:6 FOUND";
                                print "g order 13, f order 6, f^-1*g*f=g^10";
                                print "Trace profile [g,f,f^2,f^3] =",
                                      [TrS(g),TrS(ff),TrS(ff^2),TrS(ff^3)];
                                print "Fix dims [g,f,f^2,f^3] =",
                                      [FixDim(g),FixDim(ff),FixDim(ff^2),FixDim(ff^3)];
                                print "============================================================";
                            end if;

                        elif exp eq 4 then
                            ff := f^-1;
                            nf := NormalFormCount(g,ff);

                            print "  exp=4, using f^-1 to get exponent 10";
                            print "  Normal-form count for <g,f^-1> =", nf;

                            if nf eq 78 then
                                localBHits +:= 1;
                                print "============================================================";
                                print "LOCAL B=13:6 FOUND";
                                print "g order 13, f^-1 order 6, (f^-1)^-1*g*(f^-1)=g^10";
                                print "Trace profile [g,f^-1,(f^-1)^2,(f^-1)^3] =",
                                      [TrS(g),TrS(ff),TrS(ff^2),TrS(ff^3)];
                                print "Fix dims [g,f^-1,(f^-1)^2,(f^-1)^3] =",
                                      [FixDim(g),FixDim(ff),FixDim(ff^2),FixDim(ff^3)];
                                print "============================================================";
                            end if;
                        end if;
                    end if;
                end if;
            end for;
        end for;

        print "  rep summary: product hits =", profileProductHits,
              "normalising hits =", profileNormalisingHits;
    end for;
end for;

print "============================================================";
print "FINAL SUMMARY";
print "Total f^3 product hits order13 trace0 =", totalProductHits;
print "Total hits where f normalises <g> =", normalisingHits;
print "Total local B=13:6 hits =", localBHits;
print "============================================================";

if localBHits eq 0 then
    print "No local B found in this Magma Online diagnostic.";
    print "Important: this is now evidence about genuine classes, not just failed random standard pairs.";
    print "If product hits are many but normalising hits are zero, the f^3 route needs stronger local centraliser control or HPC.";
end if;

print "DONE.";
