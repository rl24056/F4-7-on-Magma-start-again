SetSeed(9);

print "============================================================";
print "TINY RANDOM PSL2(13) FINGERPRINT";
print "No Elements(). No BFS.";
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
print "Ngens(H) =", Ngens(H);

found := false;

for trial in [1..20000] do
    x := Random(H);
    y := Random(H);

    ox := Order(x);
    oy := Order(y);

    if ox mod 2 eq 0 and oy mod 3 eq 0 then
        a := x^(ox div 2);
        b := y^(oy div 3);

        if Order(a) eq 2 and Order(b) eq 3 and Order(a*b) eq 13 then
            K := sub< H | a,b >;

            if #K eq 1092 then
                print "FOUND generating PSL2(13) pair at trial", trial;
                print "Order(a), Order(b), Order(a*b) =", Order(a), Order(b), Order(a*b);
                print "Fingerprint:";
                print FP(a,b);
                found := true;
                break;
            end if;
        end if;
    end if;

    if trial mod 1000 eq 0 then
        print "trial", trial;
    end if;
end for;

print "found =", found;
print "DONE.";
