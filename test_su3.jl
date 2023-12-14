#!/usr/bin/env julia

include("su3.jl")

r1 = SU3(0, 0);  # singlet
println("SU(3) singlet")
display(r1)
@test dimension(r1) == 1
@test dual(r1) == r1

r3 = SU3(1, 0);  # quark
println("SU(3) fundamental 'quark'")
display(r3)
@test dimension(r3) == 3

r3b = SU3(1, 1);  # antiquark
println("SU(3) dual fundamental 'antiquark' / 2 fermions")
display(r3b)
@test dimension(r3b) == 3
@test dual(r3) == r3b
@test dual(r3b) == r3

r6 = SU3(2, 0)  # 2 bosons
println("SU(3) 2-bosons")
display(r6)
@test dimension(r6) == 6

r6b = SU3(2, 2)  # 2 bosons
println("SU(3) dual 2-bosons")
display(r6b)
@test dual(r6) == r6b
@test dual(r6b) == r6
@test dimension(r6b) == 6

r8 = SU3(2, 1)  # adjoint "gluon"
println("SU(3) adjoint 'gluon'")
display(r8)
@test dual(r8) == r8
@test dimension(r8) == 8


# bruteforce test dimensions in product
for r1A = 0:10
    for r2A = 0:r1A
        rA = SU3(r1A, r2A)
        dA = dimension(rA)
        for r1B = 0:10
            for r2B = 0:r1B
                rB = SU3(r1B, r2B)
                degen, irreps = rA ⊗ rB
                dim = sum(d * dimension(irr) for (d, irr) in zip(degen, irreps))
                @test dim == dA * dimension(rB)
            end
        end
    end
end


println()
# quark ⊗ antiquark = meson ⊕ gluon
degen, irreps = r3 ⊗ r3b
println("3 ⊗ 3b =")
for (d, irr) in zip(degen, irreps)
    println(d, "⋅", dimension(irr))
    display(irr)
end

println()
# fusion ring has inner degeneracies: 8 ⊗ 8 = 1 ⊕ 2⋅8 ⊕ 10 ⊕ 10b ⊕ 27
degen, irreps = r8 ⊗ r8
println("8 ⊗ 8 =")
for (d, irr) in zip(degen, irreps)
    println(d, "⋅", dimension(irr))
    display(irr)
end
