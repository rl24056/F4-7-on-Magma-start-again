SetSeed(6);

print "============================================================";
print "SHORT F4(7) SEARCH: FORCE INVOLUTION TRACE 1";
print "Search |a|=2, Tr(a)=1, |b|=3, Tr(b)=+-1, |ab|=13, Tr(ab)=0";
print "Then enumerate only if product hit appears.";
print "============================================================";

TRIALS := 120000;
REP_TRIALS := 5000;
PROGRESS := 5000;
ENUM_CAP := 1092;

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
        if x^(n div p) eq Id then return false; end if;
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

/* ---------------- build G ---------------- */

F := GF(7);
G := ChevalleyGroup("F",4,F);
Id := Identity(G);

print "Degree(G) =", Degree(G);
print "BaseRing(G) =", BaseRing(G);

/* ---------------- build small pools ---------------- */

APool := [ G | ];      // involutions trace 1 only
BPool := [ G | ];      // order 3 trace -1 or 1

for t in [1..REP_TRIALS] do
    r := Random(G);
    o := Order(r);

    if o mod 2 eq 0 then
        a := r^(o div 2);
        if IsExactOrder(a,2,Id) and TrS(a) eq 1 then
            AddIfNew(~APool,a);
        end if;
    end if;

    if o mod 3 eq 0 then
        b := r^(o div 3);
        if IsExactOrder(b,3,Id) and (TrS(b) eq -1 or TrS(b) eq 1) then
            AddIfNew(~BPool,b);
        end if;
    end if;

    if #APool ge 20 and #BPool ge 20 then
        break;
    end if;
end for;

print "#APool trace 1 involutions =", #APool;
print "#BPool order 3 =", #BPool;

if #APool eq 0 or #BPool eq 0 then
    error "Failed to build pools. Increase REP_TRIALS.";
end if;

/* ---------------- main search ---------------- */

hits := 0;
smallHits := 0;
final := false;

for trial in [1..TRIALS] do
    a0 := APool[Random(1,#APool)];
    b0 := BPool[Random(1,#BPool)];

    ca := Random(G);
    cb := Random(G);

    a := a0^ca;
    b := b0^cb;

    g := a*b;

    if g ne Id and g^13 eq Id and TrS(g) eq 0 then
        hits +:= 1;

        print "------------------------------------------------------------";
        print "PRODUCT HIT at trial", trial;
        print "Trace(a), Trace(b), Trace(g) =", TrS(a), TrS(b), TrS(g);
        print "Order(a), Order(b), Order(ab) =", Order(a), Order(b), Order(g);

        print "Enumerating <a,b> up to 1092...";
        Hlist, closed := EnumUpTo([a,b], ENUM_CAP);

        print "closed =", closed;
        print "size/enumerated =", #Hlist;

        if closed and #Hlist eq 1092 then
            smallHits +:= 1;
            print "SUCCESS: <a,b> has order 1092.";

            /* Extract f with g^f = g^10 */
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
                nf := NormalFormCount(g,f);

                if nf eq 78 then
                    print "Clean B=<g,f> found.";
                    print "Trace(f), Trace(f^2), Trace(f^3) =", TrS(f), TrS(f^2), TrS(f^3);

                    /* Search t */
                    for t in Hlist do
                        if IsExactOrder(t,2,Id) and t^-1*f*t eq f^-1 then
                            /* t not in B check by normal forms */
                            inB := false;
                            for i in [0..12] do
                                for j in [0..5] do
                                    if t eq g^i*f^j then
                                        inB := true;
                                    end if;
                                end for;
                            end for;

                            if not inB then
                                for i in [0..12] do
                                    for j in [0..5] do
                                        q := (g^i*f^j)*t;

                                        if q ne Id and q^3 eq Id then
                                            print "============================================================";
                                            print "FINAL WALTON TRIPLE FOUND";
                                            print "|H| = 1092";
                                            print "|B| = 78";
                                            print "Order(g), Order(f), Order(t) =", Order(g), Order(f), Order(t);
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
                        end if;

                        if final then break; end if;
                    end for;
                end if;

                if final then break; end if;
            end for;
        end if;
    end if;

    if final then break; end if;

    if trial mod PROGRESS eq 0 then
        print "trial", trial, "product hits", hits, "small hits", smallHits;
    end if;
end for;

print "============================================================";
print "DONE";
print "product hits =", hits;
print "small PSL2-size hits =", smallHits;
print "final Walton triple found =", final;
print "============================================================";
