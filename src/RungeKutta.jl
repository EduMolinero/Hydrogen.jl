module RungeKutta

struct RK_Butcher_tableau{T}
    # https://www.johndcook.com/blog/2020/02/13/runge-kutta-methods/
    # see wiki page for details
    A :: Array{T,2}
    c :: Array{T,1}
    b :: Array{T,1}
    bstar :: Array{T,1}
    bool_adaptative :: Bool
    s :: Int # order of the method
end

function construct_RK_Butcher_tableau(A,c,b,bstar, type_number )
    bool_adaptative = true
    if sum(abs.(bstar)) ≈ 0.0
        bool_adaptative = false
    end
    @assert size(A,1) == size(A,2) && size(A,1) == size(c,1) && size(c,1) == size(b,1) && size(bstar,1) == size(b,1)
    s = size(c,1)
    return RK_Butcher_tableau{type_number}(map(type_number,A),map(type_number,c),map(type_number,b),map(type_number,bstar),bool_adaptative , s )
end


# Construct the RK4_tableau
function construct_Euler_tableau(T::Type = Float64)
    A = zeros(T,1,1)
    c = zeros(T,1)
    b = zeros(T,1)
    A[1,1] = convert(T, 0)
    c[1] = convert(T, 0)
    b[1] = convert(T, 1)
    bstar = b.*0
    return construct_RK_Butcher_tableau(A,c,b,bstar, T )
end
Euler_tableau = construct_Euler_tableau()

# Construct the RK4_tableau
function construct_RK4_tableau(T::Type = Float64)
    A = [0 0 0 0
         1//2 0 0 0
         0 1//2 0 0
         0 0 1 0]
    c = [0; 1 // 2; 1 // 2; 1]
    b = [1 // 6; 1 // 3; 1 // 3; 1 // 6]
    bstar = b.*0
    return construct_RK_Butcher_tableau(A,c,b,bstar, T )
end
RK4_tableau = construct_RK4_tableau()

# Construct the DP5_tableau
# Check the paper and the https://github.com/SciML/DiffEqDevTools.jl/blob/master/src/ode_tableaus.jl
function construct_DP5_tableau(T::Type = Float64)
    A = zeros(T, 7, 7)
    c = zeros(T, 7)
    b = zeros(T, 7)
    bstar = zeros(T, 7)
    c[1] = convert(T, 0)
    c[2] = convert(T, 1 // 5)
    c[3] = convert(T, 3 // 10)
    c[4] = convert(T, 4 // 5)
    c[5] = convert(T, 8 // 9)
    c[6] = convert(T, 1)
    c[7] = convert(T, 1)
    A[2, 1] = convert(T, 1 // 5)
    A[3, 1] = convert(T, 3 // 40 )
    A[3, 2] = convert(T, 9 // 40    )
    A[4, 1] = convert(T, 44 // 45   )
    A[4, 2] = convert(T, -56//15  )
    A[4, 3] = convert(T, 32//9 )
    A[5, 1] = convert(T,  19372//6561)
    A[5, 2] = convert(T, -25360//2187       )
    A[5, 3] = convert(T, 64448//6561	 )
    A[5, 4] = convert(T, -212//729    )
    A[6, 1] = convert(T, 9017//3168       )
    A[6, 2] = convert(T,  -355//33      )
    A[6, 3] = convert(T,  46732//5247      )
    A[6, 4] = convert(T,  49//176      )
    A[6, 5] = convert(T,   -5103//18656     )
    A[7, 1] = convert(T, 35//384       )
    A[7, 2] = convert(T,   0     )
    A[7, 3] = convert(T,  500//1113      )
    A[7, 4] = convert(T,  125//192      )
    A[7, 5] = convert(T,  -2187//6784     )
    A[7, 6] = convert(T,  11//84     )
    b[1] = convert(T,  35//384     )
    b[2] = convert(T,  0     )
    b[3] = convert(T,  500//1113     )
    b[4] = convert(T,  125//192     )
    b[5] = convert(T,  -2187/6784     )
    b[6] = convert(T,  11//84    )
    b[7] = convert(T,  0     )
    bstar[1] = convert(T,  5179//57600     )
    bstar[2] = convert(T,  0     )
    bstar[3] = convert(T,  7571//16695	     )
    bstar[4] = convert(T,  393//640     )
    bstar[5] = convert(T,  -92097//339200     )
    bstar[6] = convert(T,  187//2100	     )
    bstar[7] = convert(T,  1//40     )
    A = map(T, A)
    b = map(T, b)
    c = map(T, c)
    bstar = map(T, bstar)
    return construct_RK_Butcher_tableau(A,c,b,bstar, Float64 )
end
DP5_tableau = construct_DP5_tableau()



# Construct the Tsit5_tableau
# Check the paper and the https://github.com/SciML/DiffEqDevTools.jl/blob/master/src/ode_tableaus.jl
# Still do not understand how this works!!!! lets try alphaeest in DiffEqDevTools to be equal to bstar
function construct_Tsit5_tableau(T::Type = Float64)
    A = zeros(T, 7, 7)
    c = zeros(T, 7)
    b = zeros(T, 7)
    bstar = zeros(T, 7)
    αEEst = zeros(T, 7)
    c[2] = convert(T, 161 // 1000)
    c[3] = convert(T, 327 // 1000)
    c[4] = convert(T, 9 // 10)
    c[5] = convert(T,
                   big".9800255409045096857298102862870245954942137979563024768854764293221195950761080302604")
    c[6] = convert(T, 1)
    c[7] = convert(T, 1)
    A[2, 1] = convert(T, 161 // 1000)
    A[3, 1] = convert(T,
                      big"-.8480655492356988544426874250230774675121177393430391537369234245294192976164141156943e-2")
    A[3, 2] = convert(T,
                      big".3354806554923569885444268742502307746751211773934303915373692342452941929761641411569")
    A[4, 1] = convert(T,
                      big"2.897153057105493432130432594192938764924887287701866490314866693455023795137503079289")
    A[4, 2] = convert(T,
                      big"-6.359448489975074843148159912383825625952700647415626703305928850207288721235210244366")
    A[4, 3] = convert(T,
                      big"4.362295432869581411017727318190886861027813359713760212991062156752264926097707165077")
    A[5, 1] = convert(T,
                      big"5.325864828439256604428877920840511317836476253097040101202360397727981648835607691791")
    A[5, 2] = convert(T,
                      big"-11.74888356406282787774717033978577296188744178259862899288666928009020615663593781589")
    A[5, 3] = convert(T,
                      big"7.495539342889836208304604784564358155658679161518186721010132816213648793440552049753")
    A[5, 4] = convert(T,
                      big"-.9249506636175524925650207933207191611349983406029535244034750452930469056411389539635e-1")
    A[6, 1] = convert(T,
                      big"5.861455442946420028659251486982647890394337666164814434818157239052507339770711679748")
    A[6, 2] = convert(T,
                      big"-12.92096931784710929170611868178335939541780751955743459166312250439928519268343184452")
    A[6, 3] = convert(T,
                      big"8.159367898576158643180400794539253485181918321135053305748355423955009222648673734986")
    A[6, 4] = convert(T,
                      big"-.7158497328140099722453054252582973869127213147363544882721139659546372402303777878835e-1")
    A[6, 5] = convert(T,
                      big"-.2826905039406838290900305721271224146717633626879770007617876201276764571291579142206e-1")
    A[7, 1] = convert(T,
                      big".9646076681806522951816731316512876333711995238157997181903319145764851595234062815396e-1")
    A[7, 2] = convert(T, 1 // 100)
    A[7, 3] = convert(T,
                      big".4798896504144995747752495322905965199130404621990332488332634944254542060153074523509")
    A[7, 4] = convert(T,
                      big"1.379008574103741893192274821856872770756462643091360525934940067397245698027561293331")
    A[7, 5] = convert(T,
                      big"-3.290069515436080679901047585711363850115683290894936158531296799594813811049925401677")
    A[7, 6] = convert(T,
                      big"2.324710524099773982415355918398765796109060233222962411944060046314465391054716027841")
    b[1] = convert(T,
                   big".9646076681806522951816731316512876333711995238157997181903319145764851595234062815396e-1")
    b[2] = convert(T, 1 // 100)
    b[3] = convert(T,
                   big".4798896504144995747752495322905965199130404621990332488332634944254542060153074523509")
    b[4] = convert(T,
                   big"1.379008574103741893192274821856872770756462643091360525934940067397245698027561293331")
    b[5] = convert(T,
                   big"-3.290069515436080679901047585711363850115683290894936158531296799594813811049925401677")
    b[6] = convert(T,
                   big"2.324710524099773982415355918398765796109060233222962411944060046314465391054716027841")
    bstar[1] = convert(T,
                       big".9468075576583945807478876255758922856117527357724631226139574065785592789071067303271e-1")
    bstar[2] = convert(T,
                       big".9183565540343253096776363936645313759813746240984095238905939532922955247253608687270e-2")
    bstar[3] = convert(T,
                       big".4877705284247615707855642599631228241516691959761363774365216240304071651579571959813")
    bstar[4] = convert(T,
                       big"1.234297566930478985655109673884237654035539930748192848315425833500484878378061439761")
    bstar[5] = convert(T,
                       big"-2.707712349983525454881109975059321670689605166938197378763992255714444407154902012702")
    bstar[6] = convert(T,
                       big"1.866628418170587035753719399566211498666255505244122593996591602841258328965767580089")
    bstar[7] = convert(T, 1 // 66)
    A = map(T, A)
    b = map(T, b)
    c = map(T, c)
    bstar = map(T, bstar)
    return construct_RK_Butcher_tableau(A,c,b,bstar, Float64 )
end
Tsit5_tableau = construct_Tsit5_tableau()

# See interpolants at https://github.com/SciML/OrdinaryDiffEq.jl/blob/571d16cc24295080348aa6303a8d49a282e2e0b4/src/dense/interpolants.jl

end