#sygnały o zadanej charakterystyce "oblicz ... dyskretnego sygnału"
          
function rozwiazanie(;
    fp::Float64 = 434.88,
    t1::Float64 = 0.92,
    N::Int = 826,
)    
g1(t) = 2*(t-floor(t+1/2)) # piłokształtny rosnący 
g2(t) = -2*(t-floor(t+1/2)) # piła opadajaca 
g3(t) = 2*(abs(2*(t+1/4) - floor((t+1/4)+1/2))) # trójkątna
g4(t) = sign(sin(2*pi*t)) # bipolarna prostokątna

x = range(start = t1, step = 1/fp, length = N)
y = [a*g(b*t + c) for t in x]

return sum(y) # suma sygnału 
return sum(y)/length(y) # srednia
return sum(y.^2)/length(y) # moc 
return sum(y.^2) # energia
return sqrt(sum(y.^2)/length(y)) # wartosc skuteczna

end

#interpolacyjnego Whittakera-Shannona

          
function rozwiazanie(;
    m::Vector{Float64} = [-1.9, -1.8987, -1.8974, -1.8961, -1.8948, -1.8935, -1.8922, -1.8909, -1.8896, -1.8883, -1.887, -1.8857, -1.8844, -1.8831, -1.8818, -1.8805, -1.8792, -1.8779, -1.8766, -1.8753, -1.874, -1.8727, -1.8714, -1.8701, -1.8688, -1.8675, -1.8662, -1.8649, -1.8636, -1.8623, -1.861, -1.8597, -1.8584, -1.8571, -1.8558, -1.8545, -1.8532, -1.8519, -1.8506, -1.8493, -1.848, -1.8467, -1.8454, -1.8441, -1.8428, -1.8415, -1.8402, -1.8389, -1.8376, -1.8363, -1.835, -1.8337, -1.8324, -1.8311, -1.8298, -1.8285, -1.8272, -1.8259, -1.8246, -1.8233, -1.822, -1.8207, -1.8194, -1.8181, -1.8168, -1.8155, -1.8142, -1.8129, -1.8116, -1.8103, -1.809, -1.8077, -1.8064, -1.8051, -1.8038, -1.8025, -1.8012, -1.7999, -1.7986, -1.7973, -1.796, -1.7947, -1.7934, -1.7921, -1.7908, -1.7895, -1.7882, -1.7869],
    s::Vector{Float64} = [0.2384, 0.2154, 0.3467, 0.3508, 0.9385, 0.9448, 0.0322, 0.5264, 0.6935, 0.3755, 0.4112, 0.8084, 0.2717, 0.436, 0.7842, 0.5032, 0.424, 0.3929, 0.5494, 0.5764, 0.2673, 0.3766, 0.2168, 0.3493, 0.6712, 0.2229, 0.9385, 0.5739, 0.4516, 0.0411, 0.2634, 0.9194, 0.6057, 0.6863, 0.3303, 0.0489, 0.9583, 0.0187, 0.8108, 0.0164, 0.0743, 0.6953, 0.9907, 0.2156, 0.7892, 0.6026, 0.7618, 0.8341, 0.1214, 0.1445, 0.0859, 0.4364, 0.9211, 0.4757, 0.5064, 0.1144, 0.4232, 0.055, 0.6362, 0.3088, 0.5513, 0.032, 0.9701, 0.2993, 0.493, 0.8706, 0.9195, 0.6015, 0.3296, 0.5251, 0.3623, 0.3513, 0.5403, 0.013, 0.0026, 0.3674, 0.2052, 0.4334, 0.2988, 0.1113, 0.088, 0.895, 0.7062, 0.7936, 0.0358, 0.3431, 0.165, 0.19],
    t::Vector{Float64} = [-1.85476, -1.88492, -1.86841, -1.86139, -1.8883, -1.83448, -1.81017, -1.89571, -1.87491, -1.86776],
)
    t_out = zeros(length(t))
    T = m[2]-m[1]
    for i in 1 : length(t)
        for n in 1:length(s)
            t_out[i] += sinc((t[i]-m[n])/T)*s[n]
        end
    end
    
    return sum(t_out)
end

#kwantyzacja x bitowa

function rozwiazanie(;
    a::Float64 = 0.016,
    b::Float64 = 0.98,
    x::Vector{Float64} = [0.94862, 0.0241, 0.07861, 0.1851, 0.94997, 0.12357, 0.37632, 0.46596, 0.81761, 0.61411, 0.53581, 0.71855, 0.21184, 0.25854, 0.11878, 0.35106, 0.41106, 0.57201, 0.3064, 0.5018, 0.78848, 0.79263, 0.28096, 0.52367, 0.20863, 0.51147, 0.33711, 0.22929, 0.9605, 0.08728, 0.75651, 0.03799, 0.96973, 0.12431, 0.39234, 0.47093, 0.89223, 0.66858, 0.24805, 0.32817, 0.11595, 0.88675, 0.53843, 0.36247, 0.61832, 0.75149, 0.67319, 0.11964, 0.65864, 0.1502, 0.87146, 0.15946, 0.08376, 0.10753, 0.29461, 0.59939, 0.40131, 0.41222, 0.70343, 0.97864, 0.94813, 0.64549, 0.6861, 0.41727, 0.15312, 0.93553, 0.93197, 0.38447, 0.79489, 0.44495, 0.70823, 0.54983, 0.31198, 0.68062, 0.75205, 0.04621, 0.08287, 0.55077, 0.31136, 0.43562, 0.50655, 0.82521, 0.89177, 0.18876, 0.10781, 0.49479, 0.01551, 0.50157, 0.39033, 0.59818],
)
    quantize(L) = x->L[argmin(abs.(x.-L))]
    N = 7 #ilosc bitow
    L = range(start = a, stop = b, length=2^N)
    q = quantize(L)
    x_q = q.(x)
    e = x - x_q
    return #... 
end


# dft 

          
function rozwiazanie7(;
    fp::Int = 987,
    x::Vector{ComplexF64} = ComplexF64[0.8 - 0.62im, 0.68 + 0.19im, -0.06 + 0.05im, -0.23 + 0.4im, 0.17 - 0.23im, 0.43 + 0.15im, -0.93 - 0.33im, 0.81 + 0.06im, -1.37 - 0.23im, -0.78 + 0.4im, 0.13 - 0.62im, 1.91 - 0.61im, -0.63 + 0.37im, -1.69 - 0.31im, 0.05 - 0.87im, 0.87 + 0.07im, -0.62 + 0.44im, -0.78 - 0.88im, 0.51 - 0.92im, 1.11 + 0.63im, 0.73 - 0.0im, -0.39 + 0.8im, -1.32 - 0.67im, 0.43 + 0.0im, 0.1 - 0.51im, 1.2 + 0.01im, -0.17 - 0.22im, 0.15 + 0.41im, -0.17 + 0.47im, -0.67 + 0.81im, -0.13 - 0.06im, -0.23 + 0.31im, -0.01 + 0.57im, 1.42 - 0.27im, 0.74 - 0.14im, -0.1 + 1.28im, 0.59 - 0.48im, -0.42 - 1.0im, 0.04 - 0.1im, 0.09 + 0.78im, 1.38 - 0.14im, -0.61 + 0.15im, -1.36 - 0.02im, 0.5 - 1.12im, 0.65 - 0.75im, 0.63 - 0.87im, 0.32 + 0.89im],
    f::Vector{Int} = [-441, -420, -63, -42, -21, 210, 462],
)
   N = length(x)
   return sum([angle.(sum([x[n+1]cispi(-2fi/fp*n) for n in 0:N-1 ])/N) for fi in f]) 
end
rozwiazanie7()

# odpowiedz impulsowa sygnału 

function rozwiazanie(;
    x::Vector{Float64} = [2.28, 2.69, 1.28, -0.29, 4.65, -3.89, -4.74, -1.47, -3.25, -1.9, -2.74, -0.2, 0.24, 4.77, 3.5, -4.63, -3.34, -3.4, -2.67, -2.58, -0.66, -3.74, -0.9, 4.46, 0.51, 1.28, 0.23, -0.87, 0.48, 2.58, -0.07, 4.61, -3.54, -1.51, 1.6, 2.35, -4.25, -1.43, -3.53, 4.81, 3.97, -4.08, -3.55, -4.12, -4.98, 1.84, 1.29, 3.4, -1.08, -2.81, -3.33, 0.1, -2.79, 0.98, -2.5, -3.68, -1.95, -0.93, 3.3, 0.99, -0.59, 4.54],
    h::Vector{Float64} = [4.07, -1.66, -1.06, 3.0, 2.17, 4.23, -0.32, -3.56, -4.19, -1.83, -3.93, -4.93, -0.65, 4.24, -1.55, 4.23, 1.41, -3.42, 4.25, -3.79, -2.32, -4.03],
)
function conv(x,h)
    N = length(x)
    M = length(h)
    K = N+M-1
    y = zeros(Float64,K)
    for i in 1:K 
        for m in 1:N
            if i-m+1 >0 && i-m+1<=M 
                y[i] += h[m]*x[n-m+1]
            end
        end
    end
    return y 
    end
    
    y = conv(x,h)
    return #...

end

# równaniem różnicowym

          
function rozwiazanie(;
    b::Vector{Float64} = [0.005886216155083775, 0.0, -0.017658648465251326, 0.0, 0.017658648465251326, 0.0, -0.005886216155083775],
    a::Vector{Float64} = [1.0, -2.8390104753828656, 4.8984221995853146, -5.0962271628207905, 3.7271905661816733, -1.637353111316897, 0.43922072794642975],
    x::Vector{Float64} = [0.65, 0.72, 0.8, 0.58, -0.31, -0.92, -0.9, 0.21, 0.87, -0.52, -0.51, -0.37, 0.66, -0.29, -0.09, 0.44, -0.94, -0.25, -0.95, 0.15, -0.74, 0.28, -0.95, -0.99, 0.33, 0.42, -0.89, 0.89, -0.03],
    L::Int = 50,
)
    K = length(a)
    M = length(b)
    N = length(x)
    y = zeros(Float64, L)

    for n in 1:L 
        for m in 1:M
            if n-m+1>0 && n-m+1 <=N 
                y[n] += b[m]*x[n-m+1]
            end
        end
        for k in 2:K
            if n-k+1 >0 && n-k+1 <=L
                y[n] -=a[k]*y[n-k+1]
            end
        end
    end
    return sum(y.^2)/length(y)
end
rozwiazanie()

#Dany jest dyskretny system liniowy niezmienny w czasie, który jest opisany Z-transmitacją H(z).
#podaj średnie przesunięcie fazowe
          
function rozwiazanie(;
    b::Vector{Float64} = [6.238698354847942e-5, 0.0, -0.00024954793419391767, 0.0, 0.0003743219012908765, 0.0, -0.00024954793419391767, 0.0, 6.238698354847942e-5],
    a::Vector{Float64} = [1.0, -3.6330270926156896, 8.463618351718203, -12.591890016179114, 14.06407276350723, -11.131225695520087, 6.6138091257372205, -2.5089273284907283, 0.6105348075612239],
    F::Vector{Float64} = [0.02, 0.05, 0.12, 0.26, 0.39],
)
    out = zeros(Float64, length(F))
        for i in 1:length(F)
            mianownik = 0
            licznik = 0 
            for m in 1:length(b)
                licznik += b[m]*cispi(2F[i])^-(m-1)
            end
            for n in 1:length(a)
                mianownik += a[n]*cispi(2F[i])^-(n-1)
            end
            out[i] = angle(licznik/mianownik)
        end
        return sum(out)/length(out)
end
rozwiazanie()
#srednie wzmocnienie 

function rozwiazanie(;
    b::Vector{Float64} = [6.238698354847942e-5, 0.0, -0.00024954793419391767, 0.0, 0.0003743219012908765, 0.0, -0.00024954793419391767, 0.0, 6.238698354847942e-5],
    a::Vector{Float64} = [1.0, -3.6330270926156896, 8.463618351718203, -12.591890016179114, 14.06407276350723, -11.131225695520087, 6.6138091257372205, -2.5089273284907283, 0.6105348075612239],
    F::Vector{Float64} = [0.02, 0.05, 0.12, 0.26, 0.39],
)

out = zeros(Float64, length(F))

    for i in 1:length(F)
        mianownik = 0
        licznik = 0
        for l in 1:length(b)
            licznik += b[l] * cispi(2*F[l])^-(l-1)
        end
        for m in 1:length(a) 
            mianownik +=a[m] * cispi(2*F[m])^-(m-1)
        end
    out[i] = abs(licznik/mianownik)
    end
    return sum(out)/length(out) 

end




# średnie przesuniecie fazowe inne dane        
function rozwiazanie(;
    zz::Vector{ComplexF64} = ComplexF64[1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im],
    pp::Vector{ComplexF64} = ComplexF64[0.8956223942668234 - 0.25133617864363955im, 0.8956223942668234 + 0.25133617864363955im, 0.802065089637981 - 0.16477104702341286im, 0.802065089637981 + 0.16477104702341286im, 0.7564436771955586 - 0.05687993784022431im, 0.7564436771955586 + 0.05687993784022431im],
    k::Float64 = 0.5777931270964739,
    F::Vector{Float64} = [0.06, 0.1, 0.21, 0.41, 0.45],
)
out = zeros(Float64, length(F))
for i in 1:length(F)
    licznik = 1
    mianownik = 1
    for l in 1:length(zz)
        licznik *= cispi(2*F[i])-zz[l]
    end
    for m in 1:length(pp)
        mianownik *= cispi(2*F[i])-pp[m]
    end
    out[i]=angle(k*licznik/mianownik)
end
return sum(out)/length(out)



end

# średnie wzmocnienie 
function rozwiazanie(;
    zz::Vector{ComplexF64} = ComplexF64[1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im],
    pp::Vector{ComplexF64} = ComplexF64[0.8956223942668234 - 0.25133617864363955im, 0.8956223942668234 + 0.25133617864363955im, 0.802065089637981 - 0.16477104702341286im, 0.802065089637981 + 0.16477104702341286im, 0.7564436771955586 - 0.05687993784022431im, 0.7564436771955586 + 0.05687993784022431im],
    k::Float64 = 0.5777931270964739,
    F::Vector{Float64} = [0.06, 0.1, 0.21, 0.41, 0.45],
)
out = zeros(Float64,length(F))
for i in 1:length(F)
    licznik = 1
    mianownik = 1
        for l in 1:length(zz)
            licznik *=cispi(2*F[i])-zz[l]
        end
        for m in 1:length(pp)
            mianownik *=cispi(2*F[i])-pp[m]
        end
    out[i]=abs(k*licznik/mianownik)
    return sum(out)/length(out) 
end






#Dany jest dyskretny system liniowy niezmienny w czasie, który jest opisany Z-transmitacją 
# zbadaj czy jest stabilny 
          
function rozwiazanie(;
    z::Vector{ComplexF64} = ComplexF64[-1.0 + 0.0im, -1.0 + 0.0im, -1.0 + 0.0im, -1.0 + 0.0im, -1.0 + 0.0im, -1.0 + 0.0im],
    p::Vector{ComplexF64} = ComplexF64[4.661130082386245 + 1.8093828658064588im, 0.9162239951193701 - 0.3556648213432879im, 0.9178658426254936 + 0.2566400603460595im, 0.9178658426254936 - 0.2566400603460595im, 0.9307294155233051 + 0.09374688066621499im, 0.9307294155233051 - 0.09374688066621499im],
    k::Float64 = 1.6348215201772598e-6,
)
    for i in p 
        if abs(i) > 0 
            return -1
        end
    end
    for i in p 
        if abs(i) == 0
            return 0
        end
    end
    return 1
end
rozwiazanie()

# filtry 

kronecker(n::Integer)::Real = ifelse(n== 0, 1, 0)

dolno(order::Integer, F0::Float64)::Vector = [2F0*sinc(2F0*n) for n in -order/2:order/2]
gorno(order::Integer, F0::Float64)::Vector = [kronecker(Int(n))-2F0*sinc(2F0*n) for n in -order/2:order/2]
pasmowoprzepust(order::Integer,F1::Float64, F2::Float64)::Vector = [2F2*sinc(2F2*n)-2F1*sinc(2F1*n) for n in -order/2:order/2]
pasmowozaporowy(order::Integer, F1::Float64, F2::Float64)::Vector = [kronecker(Int(n))-(2F2*sinc(2F2*n) - 2F1*sinc(2F1*n)) for n in -order/2:order/2]


#okna 

 triang(M::Integer)::AbstractVector{<:Real} = [1-abs(n)/(M+1) for n in -M:M]
 hanning(M::Integer)::AbstractVector{<:Real} = [0.5(1+cos(2*pi*n/(2M+1))) for n in -M:M]
 hamming(M::Integer)::AbstractVector{<:Real} = [0.54 + 0.46*cos(2*pi*n/(2M+1)) for n in -M:M]
 blackman(M::Integer)::AbstractVector{<:Real} = [0.42 + 0.5*cos(2pi*n/(2M+1))+0.08*cos(4pi*n/(2M+1)) for n in -M:M]




function rozwiazanie2(;
    order:: Int = 86,
    fp::Float64 = 192.0,
    f0::Float64 = 32.64,
    z::Vector{Int} = [26, 24, 82],

)
a = gorno(order, f0/fp)
a = a.*blackman(Int(order/2))
a_z = [a[i] for i in z]
return sum(a_z)
end
rozwiazanie()

