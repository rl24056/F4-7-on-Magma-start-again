SetSeed(8);

print "============================================================";
print "SHORT F4(7) (2,3,13) SEARCH WITH HARDCODED PSL2(13) FINGERPRINT";
print "Use trace-2 involutions only; trace-1 gave zero product hits.";
print "No abstract PSL2 code. No Normalizer.";
print "============================================================";

TRIALS := 200000;
REP_TRIALS := 5000;
PROGRESS := 5000;
ENUM_CAP := 1092;

/* PSL2(13) fingerprint for any generating pair |a|=2, |b|=3, |ab|=13 */
TARGET_FP := [
    13,13,
    7,7,
    6,6,
    6,6,
    7,7,
    7,7,
    7,7,
    6,6,
    7,7,
    6,6
];

/* ---------------- helpers ---------------- */

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

function Comm(a,b)
    return a^-1*b^-1*a*b;
end function;

function FP(a,b)
    ab := a*b;
    ab2 := a*b^2;

    w1 := a*b*a*b^2;
    w2 := a*b^2*a*b;

    c1 := Comm(a,b);
    c2 := Comm(a,b^2);

    return [
        Order(ab), Order(ab2),
        Order(w1), Order(w2),
        Order(ab^2*ab2), Order(ab2^2*ab),
        Order(w1*ab), Order(w2*ab2),
        Order(w1^2), Order(w2^2),
        Order(ab^3*ab2), Order(ab2^3*ab),
        Order(c1), Order(c2),
        Order(c1*ab), Order(c1*ab2),
        Order(c1^2*ab), Order(c1^2*ab2),
        Order(ab*ab2*ab), Order(ab2*ab*ab2)
    ];
end function;

function EnumUpTo(gens,cap)
    P := Parent(gens[1]);
    Id := Identity(P);

    S := [ P | ];
    for x in gens do
        Append(~S,x);
        Append(~S,x^-1);
    end for;

    L := [ P | Id ];
    i := 1;

    while i le #L do
        x := L[i];

        for s in S do
            y := x*s;

            if not InList(y,L) then
                Append(~L,y);

                if #L gt cap then
                    return L,false;
                end if;
            end if;
        end for;

        i +:= 1;
    end while;

    return L,true;
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

function InB13by6(x,g,f)
    for i in [0..12] do
        for j in [0..5] do
            if x eq g^i*f^j then
                return true;
            end if;
        end for;
    end for;

    return false;
end function;

/* ---------------- build genuine F4(7) ---------------- */

F := GF(7);
G := ChevalleyGroup("F",4,F);
Id := Identity(G);

print "Degree(G) =", Degree(G);
print "BaseRing(G) =", BaseRing(G);

/* ---------------- build pools ---------------- */

APool := [ G | ];  // trace-2 involutions only
BPool := [ G | ];  // order-3 elements, both trace classes

for t in [1..REP_TRIALS] do
    r := Random(G);
    o := Order(r);

    if o mod 2 eq 0 then
        a := r^(o div 2);
        if IsExactOrder(a,2,Id) and TrS(a) eq 2 then
            AddIfNew(~APool,a);
        end if;
    end if;

    if o mod 3 eq 0 then
        b := r^(o div 3);
        if IsExactOrder(b,3,Id) and (TrS(b) eq -1 or TrS(b) eq 1) then
            AddIfNew(~BPool,b);
        end if;
    end if;

    if #APool ge 25 and #BPool ge 40 then
        break;
    end if;
end for;

print "#APool trace-2 involutions =", #APool;
print "#BPool order-3 elements =", #BPool;

if #APool eq 0 or #BPool eq 0 then
    error "Pool construction failed.";
end if;

/* ---------------- main search ---------------- */

productHits := 0;
fingerprintHits := 0;
smallHits := 0;
final := false;

HitTraceCounts := AssociativeArray();

for trial in [1..TRIALS] do
    a0 := APool[Random(1,#APool)];
    b0 := BPool[Random(1,#BPool)];

    a := a0^Random(G);
    b := b0^Random(G);

    g := a*b;

    if g ne Id and g^13 eq Id and TrS(g) eq 0 then
        productHits +:= 1;

        trkey := Sprint(<TrS(a),TrS(b)>);
        if IsDefined(HitTraceCounts,trkey) then
            HitTraceCounts[trkey] +:= 1;
        else
            HitTraceCounts[trkey] := 1;
        end if;

        fp := FP(a,b);

        if fp eq TARGET_FP then
            fingerprintHits +:= 1;

            print "------------------------------------------------------------";
            print "PSL2 fingerprint hit at trial", trial;
            print "productHits =", productHits;
            print "fingerprintHits =", fingerprintHits;
            print "Trace(a),Trace(b),Trace(g) =", TrS(a),TrS(b),TrS(g);
            print "fingerprint =", fp;

            print "Enumerating <a,b> up to 1092...";
            Hlist, closed := EnumUpTo([a,b], ENUM_CAP);

            print "closed =", closed;
            print "size/enumerated =", #Hlist;

            if closed and #Hlist eq 1092 then
                smallHits +:= 1;
                print "SUCCESS: <a,b> has order 1092.";

                /* extract f */
                FCands := [ G | ];

                for x in Hlist do
                    if IsExactOrder(x,6,Id) then
                        e := ActionExponent(g,x);

                        if e eq 10 then
                            Append(~FCands,x);
                        elif e eq 4 then
                            Append(~FCands,x^-1);
                        end if;
                    end if;
                end for;

                print "#f candidates =", #FCands;

                for f in FCands do
                    if NormalFormCount(g,f) eq 78 then
                        print "Clean B=<g,f> found.";
                        print "Trace(f),Trace(f^2),Trace(f^3) =", TrS(f),TrS(f^2),TrS(f^3);

                        for t in Hlist do
                            if IsExactOrder(t,2,Id) and t^-1*f*t eq f^-1 and not InB13by6(t,g,f) then
                                for i in [0..12] do
                                    for j in [0..5] do
                                        q := (g^i*f^j)*t;

                                        if q ne Id and q^3 eq Id then
                                            print "============================================================";
                                            print "FINAL WALTON TRIPLE FOUND";
                                            print "|H| = 1092";
                                            print "|B| = 78";
                                            print "Order(g),Order(f),Order(t) =", Order(g),Order(f),Order(t);
                                            print "Trace profile [g,f,f^2,f^3,t,g*t] =",
                                                  [TrS(g),TrS(f),TrS(f^2),TrS(f^3),TrS(t),TrS(g*t)];
                                            print "f^-1*g*f = g^10 ?", f^-1*g*f eq g^10;
                                            print "t^-1*f*t = f^-1 ?", t^-1*f*t eq f^-1;
                                            print "Bruhat witness i,j =", i,j;
                                            print "Order((g^i f^j)t) =", Order(q);
                                            print "============================================================";

                                            final := true;
                                            break;
                                        end if;
                                    end for;

                                    if final then break; end if;
                                end for;
                            end if;

                            if final then break; end if;
                        end for;
                    end if;

                    if final then break; end if;
                end for;
            end if;
        end if;
    end if;

    if final then break; end if;

    if trial mod PROGRESS eq 0 then
        print "trial", trial;
        print "productHits =", productHits;
        print "fingerprintHits =", fingerprintHits;
        print "smallHits =", smallHits;

        print "product-hit trace counts:";
        for k in Keys(HitTraceCounts) do
            print k, HitTraceCounts[k];
        end for;
    end if;
end for;

print "============================================================";
print "DONE";
print "productHits =", productHits;
print "fingerprintHits =", fingerprintHits;
print "small PSL2-size hits =", smallHits;
print "final Walton triple found =", final;
print "============================================================";

if fingerprintHits eq 0 then
    print "No candidate matched the PSL2(13) word fingerprint.";
    print "This strongly suggests the random standard-pair route is too sparse for Magma Online.";
elif smallHits eq 0 then
    print "Some candidates matched the PSL2 word fingerprint but still generated large groups.";
    print "Then we need stronger relations or HPC.";
end if;
