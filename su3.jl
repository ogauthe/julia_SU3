using Printf
using Test

"""
Representation theory for the Lie group SU(3), aka A_2

Irreducible representations, or irreps, are labelled their Young tableau [row1, row2], with
row1 and row2 integer such that row1 >= row2 >= 0. The highest weight (row1 - row2, row2)
is also used in the litterature.

An irrep [row1, row2] can be understood as the many-body Hilbert space containing two
different species of bosons, with row1 - row2 bosonic quarks and row2 bosonic antiquarks.
One antiquark can be seen as the product of 2-fermionic quarks.
"""
struct SU3
    row1::Int
    row2::Int
    SU3(n1, n2) = n1 >= n2 >= 0 ? new(n1, n2) : error("Invalid Young tableau")
end

"""
    highest_weight(yt)

Compute the irrep unique highest weight.
"""
highest_weight(yt::SU3) = (yt.row1 - yt.row2, yt.row2)

"""
    nboxes(yt)

Compute the number of boxes in the Young tableau.
"""
nboxes(yt::SU3) = yt.row1 + yt.row2

# SU(3) has a shortcut for hook length formula
dimension(yt::SU3) = (yt.row1 - yt.row2 + 1) * (yt.row2 + 1) * (yt.row1 + 2) ÷ 2

Base.show(io::IO, yt::SU3) = @printf(io, "[%d, %d]", yt.row1, yt.row2)

# display pretty Young tableau with utf8 box char
function Base.show(io::IO, ::MIME"text/plain", yt::SU3)
    if yt.row1 == 0
        println("●")
        return
    end
    println("┌─" * "┬─" ^ (yt.row1 - 1) * "┐")
    if yt.row2 == 0
        println("└─" * "┴─" ^ (yt.row1 - 1) * "┘")
        return
    end
        println(
            "├─",
            "┼─" ^ (yt.row2 - 1 + (yt.row1 > yt.row2)),
            "┴─" ^ max(0, (yt.row1 - yt.row2 - 1)),
            "┤" ^ (yt.row2 == yt.row1),
            "┘" ^ (yt.row1 > yt.row2)
        )
        println("└─" * "┴─"^ (yt.row2 - 1) * "┘")
    return
end

"""
    dual(yt)

Return the dual irrep, which is usually not equivalent. This operation corresponds to
swapping quarks and antiquarks. Sometimes called "charge conjugation" or "electron-hole
symmetry"
"""
dual(yt::SU3) = SU3(yt.row1, yt.row1 - yt.row2)

"""
Compute SU(3) fusion rules using Littlewood-Richardson rule for Young tableaus.
See e.g. Di Francesco, Mathieu and Sénéchal, section 13.5.3.
"""
function ⊗(left::SU3, right::SU3)
    if nboxes(right) > nboxes(left)
        return right ⊗ left
    end
    if right.row1 == 0
        return [1], [left]
    end

    out = []

    # put a23 boxes on 2nd or 3rd line
    a23max1 = 2 * left.row1  # row2a <= row1a
    a23max2 = right.row1  # a2 + a3 <= total number of a
    a23max = min(a23max1, a23max2)
    for a23 = 0:a23max
        a3min1 = left.row2 + 2 * a23 - left.row1 - right.row1
        a3min2 = left.row2 - left.row1 + a23  # no a below a: row2a <= row1
        a3min = max(0, a3min1, a3min2)
        a3max1 = left.row2  # row3a <= row2a
        a3max2 = a23  # a3 <= a2 + a3
        a3max3 = right.row1 - right.row2  # more a than b, right to left: b2 + b3 <= a1 + a2
        a3max = min(a3max1, a3max2, a3max3)
        for a3 = a3min: a3max
            a2 = a23 - a3
            row1a = left.row1 + right.row1 - a23
            row2a = left.row2 + a23 - a3

            # cannot put any b on 1st line: row1ab = row1a
            b3min1 = row2a + right.row2 - row1a  # row2ab <= row1ab = row1a
            b3min2 = right.row2 + a23 - right.row1
            b3min = max(0, b3min1, b3min2)
            b3max1 = right.row2  # only other.row2 b boxes
            b3max2 = (row2a + right.row2 - a3) // 2  # row3ab >= row2ab
            b3max3 = right.row1 - a3  # more a than b, right to left: b2 <= a1
            b3max4 = row2a - a3  # no b below b: row2a >= row3ab
            b3max = min(b3max1, b3max2, b3max3, b3max4)
            for b3 = b3min: b3max
                b2 = right.row2 - b3
                row2ab = row2a + b2
                row3ab = a3 + b3
                yt = SU3(row1a - row3ab, row2ab - row3ab)

                push!(out, yt)
            end
        end
    end

    irreps = sort(unique(out), by=yt->(dimension(yt), nboxes(yt), yt.row2))
    degen = [count(==(irr), out) for irr in irreps]

    return degen, irreps
end
