SetSeed(10);

print "============================================================";
print "PSL2(13) MULTI-FINGERPRINT CHECK";
print "Check whether generating (2,3,13) pairs have one or many fingerprints.";
print "============================================================";

function FP(a,b)
    ab := a*b;
    ab2 := a*b^2;
    w1 := a*b*a*b^2;
    w2 := a*b^2*a*b;
    c1 := a^-1*b^-1*a*b;
    c2 := a^-1*(b^2)^-1*a*b^2;

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

H := PSL(2,13);
print "#H =", #H;

fps := AssociativeArray();
pairs := 0;
gens := 0;

for trial in [1..5000] do
    x := Random(H);
    y := Random(H);

    ox := Order(x);
    oy := Order(y);

    if ox mod 2 eq 0 and oy mod 3 eq 0 then
        a := x^(ox div 2);
        b := y^(oy div 3);

        if Order(a) eq 2 and Order(b) eq 3 and Order(a*b) eq 13 then
            pairs +:= 1;

            K := sub< H | a,b >;

            if #K eq 1092 then
                gens +:= 1;
                key := Sprint(FP(a,b));

                if IsDefined(fps,key) then
                    fps[key] +:= 1;
                else
                    fps[key] := 1;
                end if;
            end if;
        end if;
    end if;

    if trial mod 500 eq 0 then
        print "trial", trial, "pairs", pairs, "generating", gens, "distinct fps", #Keys(fps);
    end if;
end for;

print "============================================================";
print "SUMMARY";
print "pairs with |a*b|=13 =", pairs;
print "generating pairs =", gens;
print "distinct fingerprints =", #Keys(fps);
print "============================================================";

for k in Keys(fps) do
    print k, fps[k];
end for;

print "DONE.";
