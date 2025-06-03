using Test
using SecondQuantizedAlgebra
using SymbolicUtils
using Symbolics

const sqa = SecondQuantizedAlgebra

@testset "average_sums" begin

N = 2
ha = NLevelSpace(Symbol(:atom),2)
hf = FockSpace(:cavity)
h = hf⊗ha

ind(i) = Index(h,i,N,ha)

g(k) = IndexedVariable(:g,k)
Γij = DoubleIndexedVariable(:Γ,ind(:i),ind(:j))
σ(i,j,k) = IndexedOperator(Transition(h,:σ,i,j),k)

σn(i,j,k) = NumberedOperator(Transition(h,:σ,i,j),k)
Ω(i,j) = IndexedVariable(:Ω,i,j;identical=false)

@test Ω(ind(:i),ind(:i)) == 0

a = Destroy(h,:a)

@test(isequal(average(2*σ(1,2,ind(:k))),2*average(σ(1,2,ind(:k)))))
@test(isequal(average(g(ind(:k))*σ(2,2,ind(:k))),g(ind(:k))*average(σ(2,2,ind(:k)))))
@test(isequal(average(g(ind(:k))),g(ind(:k))))

sum1 = SingleSum(σ(1,2,ind(:k)),ind(:k))
σn(i,j,k) = NumberedOperator(Transition(h,:σ,i,j),k)
@test(isequal(σn(1,2,1)+σn(2,1,1),NumberedOperator(Transition(h,:σ,1,2)+Transition(h,:σ,2,1),1)))

#test insert_index
@test(isequal(σn(2,2,1),insert_index(σ(2,2,ind(:j)),ind(:j),1)))
@test(isequal(σ(1,2,ind(:j)),insert_index(σ(1,2,ind(:j)),ind(:k),2)))
@test(isequal(1,insert_index(1,ind(:k),1)))

sum2 = average(sum1*σ(1,2,ind(:l)))

@test(!isequal(σn(2,2,1),insert_index(sum2,ind(:j),1)))


gamma = insert_index(Γij,ind(:i),1)
gamma2 = insert_index(Γij,ind(:j),2)
gamma2_ = insert_index(Γij,ind(:i),1)
g_ = insert_index(g(ind(:j)),ind(:j),1)
@test g_ isa SymbolicUtils.BasicSymbolic
@test gamma isa SymbolicUtils.BasicSymbolic{sqa.DoubleNumberedVariable}

gamma_ = insert_index(gamma,ind(:j),2)
# @test gamma_ isa SymbolicUtils.BasicSymbolic{Parameter}
@test gamma_ isa SymbolicUtils.BasicSymbolic{Complex{Real}}

@test !isequal(gamma,gamma_)
@test !isequal(gamma,g_)
@test isequal(gamma2_,gamma)
@test !isequal(gamma,gamma2)

sumterm = σ(1,2,ind(:i))*σ(2,1,ind(:j))*σ(2,2,ind(:k))
inds = sqa.get_indices(sumterm)
@test isequal([ind(:i),ind(:j),ind(:k)],inds)

pind = Index(h,:p,5,ha)
@test isequal(4*a,Σ(a,pind,[ind(:i)]))
@test isequal(average(Σ(σ(1,2,ind(:i)),ind(:i))),sqa.IndexedAverageSum(average(σ(1,2,ind(:i))),ind(:i),[]))

avrgTerm = average(Σ(σ(2,1,ind(:i))*σ(1,2,ind(:j)),ind(:i)))
@test avrgTerm isa SymbolicUtils.BasicSymbolic && operation(avrgTerm) === +
# ADsum1 = sqa.IndexedAverageDoubleSum(avrgTerm,ind(:j),[ind(:i)])
# @test ADsum1 isa SymbolicUtils.BasicSymbolic && operation(ADsum1) === +
# @test SymbolicUtils.metadata(arguments(ADsum1)[it_2])[sqa.IndexedAverageDoubleSum] isa sqa.IndexedAverageDoubleSum
# @test SymbolicUtils.metadata(arguments(ADsum1)[2])[sqa.IndexedAverageSum] isa sqa.IndexedAverageSum

# executing it line by line gives a different order of the arguments?? arguments(ADsum1)
ADsum1 = simplify(sqa.IndexedAverageDoubleSum(avrgTerm,ind(:j),[ind(:i)]))
it_1 = findfirst(x->typeof((x))==SymbolicUtils.BasicSymbolic{IndexedAverageSum}, arguments(ADsum1))
it_2 = findfirst(x->typeof((x))==SymbolicUtils.BasicSymbolic{IndexedAverageDoubleSum}, arguments(ADsum1))
@test ADsum1 isa SymbolicUtils.BasicSymbolic && operation(ADsum1) === +
@test SymbolicUtils.metadata(arguments(ADsum1)[it_2])[sqa.IndexedAverageDoubleSum] isa sqa.IndexedAverageDoubleSum
@test SymbolicUtils.metadata(arguments(ADsum1)[it_1])[sqa.IndexedAverageSum] isa sqa.IndexedAverageSum


@test isequal(sqa.SpecialIndexedAverage(average(σ(1,2,ind(:i))),[(ind(:i),ind(:j))])+sqa.SpecialIndexedAverage(average(σ(2,1,ind(:j))),[(ind(:i),ind(:j))]),
sqa.SpecialIndexedAverage(average(σ(1,2,ind(:i))) + average(σ(2,1,ind(:j))),[(ind(:i),ind(:j))]))

@test sqa.SpecialIndexedAverage(average(0),[(ind(:i),ind(:j))]) == 0
@test SymbolicUtils.metadata(sqa.SpecialIndexedAverage(average(σ(2,1,ind(:i))),[(ind(:i),ind(:j))]))[sqa.SpecialIndexedAverage] isa sqa.SpecialIndexedAverage

@test sqa.undo_average(arguments(ADsum1)[it_2]) isa sqa.DoubleSum
@test isequal(simplify(Σ(Σ(σ(2,1,ind(:i))*σ(1,2,ind(:j)),ind(:i)),ind(:j),[ind(:i)])),simplify(sqa.undo_average(ADsum1)))

@test σ(1,2,ind(:i))*σ(2,1,ind(:j))*σn(2,2,3) isa sqa.QMul
@test σn(2,2,3)*σ(1,2,ind(:i))*σ(2,1,ind(:j)) isa sqa.QMul

# @test SymbolicUtils.iscall(sum_A.metadata) == false
@test sqa.IndexedAverageSum(1) == 1

specAvrg = sqa.SpecialIndexedAverage(average(σ(2,1,ind(:i))*σ(1,2,ind(:j))),[(ind(:i),ind(:j))])

@test isequal("(i≠1)",sqa.SecondQuantizedAlgebra.writeNeqs([(ind(:i),1)]))
@test (isequal(SymbolicUtils.arguments(SymbolicUtils.arguments(ADsum1)[it_2]),SymbolicUtils.arguments(avrgTerm)[1]) || isequal(SymbolicUtils.arguments(SymbolicUtils.arguments(ADsum1)[it_2]),SymbolicUtils.arguments(avrgTerm)[2]) )
# SymbolicUtilsv1.4.0 argument order changed
@test isequal(SymbolicUtils.arguments(SymbolicUtils.arguments(SymbolicUtils.arguments(ADsum1)[it_2])),SymbolicUtils.arguments(average(σ(2,1,ind(:i))*σ(1,2,ind(:j)))))
@test isequal(SymbolicUtils.arguments(specAvrg),SymbolicUtils.arguments(average(σ(2,1,ind(:i))*σ(1,2,ind(:j)))))

@test isequal(sqa.insert_index(σ(1,2,ind(:j))*σn(1,2,2),ind(:j),1),sqa.insert_index(sqa.insert_index(σ(1,2,ind(:i))*σ(1,2,ind(:j)),ind(:i),2),ind(:j),1))

@test isequal(sqa.IndexedAverageSum(g(ind(:i)),ind(:i),[]),average(Σ(g(ind(:i)),ind(:i),[])))
@test sqa.IndexedAverageSum(g(ind(:i)),ind(:i),[]) isa SymbolicUtils.BasicSymbolic{sqa.IndexedAverageSum}

@test ind(:i) ∈ sqa.get_indices(g(ind(:i)))

@cnumbers N_
ind2(i) = Index(h,i,N_,ha)

#functions for checking if indices occure in specific terms
function containsIndexedOps(term::Average)
    arg_ = arguments(term)
    if arg_[1] isa sqa.QMul
        for arg in arg_[1].args_nc
            if arg isa sqa.IndexedOperator
                return true
            end
        end
    else
        return arg_[1] isa sqa.IndexedOperator
    end
    return false
end
containsIndex(term::Average,ind::Index) = ind ∈ get_indices(term)

@test containsIndexedOps(average(a*σ(2,1,ind(:i))*σ(1,2,ind(:j))))
@test !(containsIndexedOps(average(a'*a)))

end
