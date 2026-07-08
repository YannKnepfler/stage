############################## AnderSim v.1 ####################################
# Incompressível - Isotérmico - Sem efeitos difusivos - Escrito todo em série ##
############################ By Anderson Viana #################################
#
include("./TDMA.jl")
using Printf
using Plots
using CSV
using DataFrames
using LinearAlgebra
#
# Leitura do Arquivo de Entrada
n_dados = 16 # número de dados de entrada
dados = zeros(n_dados) # armazenador dos dados de entrada
texto = readlines("Entrada.txt")
for i = 2:2:2*n_dados
    aux = texto[i]
    aux2 = Int(i / 2)
    global dados[aux2] = parse(Float64, aux)
end
#
# Variáveis de Entrada
Comp = dados[1]
D = dados[2]
θ = dados[3]
n_segmentos = dados[4]
tempo = dados[5]
uL_inlet = dados[6]
uG_inlet = dados[7]
ρL = dados[8]
ρG = dados[9]
P_out = dados[10]
αG_inlet = dados[11]
αG_i = dados[12]
CFL = dados[13]
Tol_P = dados[14]
Tol_u = dados[15]
Tol_α = dados[16]
#
# Parâmetros
N = Int(n_segmentos) # Número de nós da malha principal
A = (pi * D^2) / 4  # Área da seção (m²) - nesta versão é constante!
αL_inlet = 1 - αG_inlet # fração de líquido na entrada
αL_i = 1 - αG_i # fração de líquido inicial
u_inlet = [uL_inlet, uG_inlet]
Δx = Comp / n_segmentos
aux = maximum(u_inlet)
Δt = (CFL * Δx) / aux
Nt = Int(tempo ÷ Δt)
#
# Variáveis
P = zeros(N)
αL = zeros(N)
αG = zeros(N)
uL = zeros(N+1)
uG = zeros(N+1)
P0 = zeros(N)
αL0 = zeros(N)
αG0 = zeros(N)
uL0 = zeros(N+1)
uG0 = zeros(N+1)
PØ = zeros(N)
αLØ = zeros(N)
αGØ = zeros(N)
uLØ = zeros(N+1)
uGØ = zeros(N+1)
#
# Condições Iniciais
uL[:] .= uL_inlet
uG[:] .= uG_inlet
αL[:] .= αL_i
αG[:] .= αG_i
P[:] .= P_out
#
# Mapa de Valores
# Iniciais
uL0[:] .= uL_inlet
uG0[:] .= uG_inlet
αL0[:] .= αL_i
αG0[:] .= αG_i
P0[:] .= P_out
# Supostos
uLØ[:] .= uL_inlet
uGØ[:] .= uG_inlet
αLØ[:] .= αL_i
αGØ[:] .= αG_i
PØ[:] .= P_out
# Corretivos
uLδ = zeros(N+1)
uGδ = zeros(N+1)
Pδ = zeros(N)
#
# Vetores para resolução dos sistemas Lineares - Velocidades
d_L = zeros(N+1) # diagonal da matriz A do sistema Ax = b    Líquido
du_L = zeros(N+1) # diagonal inferior
dl_L = zeros(N+1) # diagonal superior
b_L = zeros(N+1) # RHS do sistema
d_G = zeros(N+1) # diagonal da matriz A do sistema Ax = b    Gás
du_G = zeros(N+1)
dl_G = zeros(N+1)
b_G = zeros(N+1)
#
# Condições de Contorno - Velocidades
d_L[1] = 1.0
du_L[1] = 0.0
dl_L[1] = 0.0
b_L[1] = uL_inlet
d_G[1] = 1.0
du_G[1] = 0.0
dl_G[1] = 0.0
b_G[1] = uG_inlet
d_L[N+1] = 1.0
du_L[N+1] = 0.0
dl_L[N+1] = -1.0
b_L[N+1] = 0 # Gradiente nulo
d_G[N+1] = 1.0
du_G[N+1] = 0.0
dl_G[N+1] = -1.0
b_G[N+1] = 0 # Gradiente nulo
#
# Vetores para resolução dos sistemas Lineares - Frações
d_Lα = zeros(N) # diagonal da matriz A do sistema Ax = b    Líquido
du_Lα = zeros(N) # diagonal inferior
dl_Lα = zeros(N) # diagonal superior
b_Lα = zeros(N) # RHS do sistema
d_Gα = zeros(N) # diagonal da matriz A do sistema Ax = b    Gás
du_Gα = zeros(N)
dl_Gα = zeros(N)
b_Gα = zeros(N)
# Condições de Contorno - Frações
d_Lα[1] = 1.0
du_Lα[1] = 0.0
dl_Lα[1] = 0.0
b_Lα[1] = αL_inlet
d_Gα[1] = 1.0
du_Gα[1] = 0.0
dl_Gα[1] = 0.0
b_Gα[1] = αG_inlet
d_Lα[N] = 1.0
du_Lα[N] = 0.0
dl_Lα[N] = -1.0
b_Lα[N] = 0 # Gradiente nulo
d_Gα[N] = 1.0
du_Gα[N] = 0.0
dl_Gα[N] = -1.0
b_Gα[N] = 0.0 # Gradiente nulo
#
# Vetores para resolução dos sistemas Lineares - Pressão
d = zeros(N) # diagonal da matriz A do sistema Ax = b
du = zeros(N) # diagonal inferior
dl = zeros(N) # diagonal superior
b = zeros(N)
# Condições de Contorno - Pressão
d[1] = 1.0
du[1] = -1.0
dl[1] = 0.0
b[1] = 0.0 # Gradiente nulo
d[N] = 1.0
du[N] = 0.0
dl[N] = 0.0
b[N] = 0.0 # Condição de Dirichlet - Não há variação - Lembrando que esse
           # sistema resolve a variação da pressão!
#
# Processamento
for t = 1:Nt
    print(t*Δt)
    print("\n")
    #
    global uL0
    global uG0
    global αL0
    global αG0
    global P0
    global uLØ
    global uGØ
    global αLØ
    global αGØ
    global PØ
    global uL
    global uG
    global αL
    global αG
    global P
    #
    # Empíricas
    τ_L = 0.0 # implementar!!!
    τ_G = 0.0 # implementar!!!
    Γ_G = 0.0 # implementar!!!
    Γ_L = -Γ_G
    #
    cont = 0
    #
    # Erros
    e_uL = 1
    e_uG = 1
    e_αL = 1
    e_αG = 1
    e_P = 1
    #
    aux_uLØ = zeros(N+1)
    aux_uGØ = zeros(N+1)
    #
    while e_uL>Tol_u || e_uG>Tol_u || e_αL>Tol_α || e_αG>Tol_α || e_P>Tol_P
        #
        # Velocidades Supostas -------------------------------------------------
        #
        e_picard_L = 1
        e_picard_G = 1
        #
        cont = 1
        while (e_picard_L > Tol_u || e_picard_G > Tol_u) # && cont < 1000
            cont = cont + 1
            #
            # Condições de Contorno de Final de Domínio
            d_L[N+1] = 1.0
            du_L[N+1] = 0.0
            dl_L[N+1] = -1.0
            b_L[N+1] = 0 # Gradiente nulo
            d_G[N+1] = 1.0
            du_G[N+1] = 0.0
            dl_G[N+1] = -1.0
            b_G[N+1] = 0 # Gradiente nulo
            #
            for i=1:(N-1)
                # Velocidade Relativa
                if Γ_G >= 0
                    Ψ = 0
                else
                    Ψ = 1
                end
                ui = Ψ * uGØ[i+1] + (1 - Ψ) * uLØ[i+1]
                #
                # Coeficientes Pós Agrupamento (Tudo em upwind)
                mt_L = (ρL * αL0[i] + ρL * αL0[i+1]) / (2 * Δt)
                mt_G = (ρG * αG0[i] + ρG * αG0[i+1]) / (2 * Δt)
                mx_L = (αLØ[i] + αLØ[i+1]) / (2 * Δx)
                mx_G = (αGØ[i] + αGØ[i+1]) / (2 * Δx)
                mg_L = 9.81*(ρL*αLØ[i] + ρL*αLØ[i+1])*sin(θ)/2
                mg_G = 9.81*(ρG*αGØ[i] + ρG*αGØ[i+1])*sin(θ)/2
                G_Les = ρL * αLØ[i+1] * (uLØ[i+1] + uLØ[i+2])/2
                G_Lws = ρL * αLØ[i] * (uLØ[i] + uLØ[i+1])/2
                G_Ges = ρG * αGØ[i+1] * (uGØ[i+1] + uGØ[i+2])/2
                G_Gws = ρG * αGØ[i] * (uGØ[i] + uGØ[i+1])/2
                # Montagem dos sitemas:
                d_L[i+1] = mt_L + G_Lws / Δx + Γ_L
                du_L[i+1] = 0.0
                dl_L[i+1] = -G_Lws / Δx
                b_L[i+1]=mt_L*uL0[i+1]-mx_L*(PØ[i+1]-PØ[i])-mg_L-Γ_L+Γ_L*ui
                #
                if uGØ[i+1] >= 0
                    d_G[i+1] = mt_G + (G_Gws/Δx) + Γ_G
                    du_G[i+1] = 0.0
                    dl_G[i+1] = -G_Gws/Δx
                else
                    d_G[i+1] = mt_G + Γ_G - G_Ges/Δx
                    du_G[i+1] = G_Ges/Δx
                    dl_G[i+1] = 0.0
                end
                b_G[i+1]=mt_G*uG0[i+1]-mx_G*(PØ[i+1]-PØ[i])-mg_G-Γ_G+Γ_G*ui

            end
            #
            aux_uLØ[:] = uLØ[:] # pra ajudar a calcular os erros
            aux_uGØ[:] = uGØ[:] # pra ajudar a calcular os erros
            #
            #
            uLØ = TDMA(d_L, du_L, dl_L, b_L)
            uGØ = TDMA(d_G, du_G, dl_G, b_G)
            #
            # uLØ[N+1] = 2uLØ[N] - uLØ[N-1]
            # uGØ[N+1] = 2*uGØ[N] - uGØ[N-1]
            #
            #
            e_picard_L = norm(uLØ-aux_uLØ)/norm(aux_uLØ)
            e_picard_G = norm(uGØ-aux_uGØ)/norm(aux_uGØ)
            #
            #
        end
        #
        # Pressão --------------------------------------------------------------
        #
        # Condições no final do domínio
        d[N] = 1.0
        du[N] = 0.0
        dl[N] = 0.0
        b[N] = 0.0
        #
        for i = 2:(N-1)
            # Coeficientes Pós Agrupamento (Tudo em upwind)
            d_Le =(Δt / Δx) * (αLØ[i] + αLØ[i+1]) / (ρL * αLØ[i] + ρL * αLØ[i+1])
            d_Lw =(Δt / Δx) * (αLØ[i] + αLØ[i-1]) / (ρL * αLØ[i] + ρL * αLØ[i-1])
            d_Ge =(Δt / Δx) * (αGØ[i] + αGØ[i+1]) / (ρG * αGØ[i] + ρG * αGØ[i+1])
            d_Gw =(Δt / Δx) * (αGØ[i] + αGØ[i-1]) / (ρG * αGØ[i] + ρG * αGØ[i-1])
            G_Le = αLØ[i] * uLØ[i+1]
            G_Lw = αLØ[i-1] * uLØ[i]
            if uG[i] >= 0
                G_Ge = αGØ[i] * uGØ[i+1]
                G_Gw = αGØ[i-1] * uGØ[i]
                alfaG_e = αGØ[i]
                alfaG_w = αGØ[i-1]
            else
                G_Ge = αGØ[i+1] * uGØ[i+1]
                G_Gw = αGØ[i] * uGØ[i]
                alfaG_e = αGØ[i+1]
                alfaG_w = αGØ[i]
            end
            # Montagem do sitema:
            d[i] = (alfaG_e * d_Ge + alfaG_w * d_Gw + αLØ[i] * d_Le +
                αLØ[i-1] * d_Lw)
            du[i] = -(alfaG_e * d_Ge + αLØ[i] * d_Le)
            dl[i] = -alfaG_w * d_Gw - αLØ[i-1] * d_Lw
            b[i] = (-(Δx * αGØ[i]) / (Δt) + (Δx * αG0[i]) / (Δt) -
                (Δx * αLØ[i]) / (Δt) + (Δx * αL0[i]) / (Δt) - G_Ge + G_Gw
                - G_Le + G_Lw)
        end
        #
        Pδ = TDMA(d, du, dl, b)
        #
        Pδ[1] = 2*Pδ[2] - Pδ[3]
        #
        # Correções --------------------------------------------------------
        for i = 1:(N-1)
            de_L =(Δt/Δx)*(αLØ[i] + αLØ[i+1])/(ρL*αLØ[i] + ρL*αLØ[i+1])
            de_G =(Δt/Δx)*(αGØ[i] + αGØ[i+1])/(ρG*αGØ[i] + ρG*αGØ[i+1])
            uLδ[i+1] = de_L * (Pδ[i] - Pδ[i+1])
            uGδ[i+1] = de_G * (Pδ[i] - Pδ[i+1])
        end
        P[:] = PØ[:] + Pδ[:] #0.3*Pδ[:]
        uL[:] = uLØ[:] + uLδ[:]  #0.7*uLδ
        uG[:] = uGØ[:] + uGδ[:]  #0.7*uGδ
        # Contorno:
        uL[N+1] = uL[N]
        uG[N+1] = uG[N]
        #
        # Alfas ------------------------------------------------------------
        #
        # condições no final do domínio
        d_Lα[N] = 1.0
        du_Lα[N] = 0.0
        dl_Lα[N] = -1.0
        b_Lα[N] = 0 # Gradiente nulo
        d_Gα[N] = 1.0
        du_Gα[N] = 0.0
        dl_Gα[N] = -1.0
        b_Gα[N] = 0.0 # Gradiente nulo
        # AGRUPAMENTO:
        for i = 2:(N-1)
            # Coeficientes Pós Agrupamento (Tudo em upwind)
            d_Lα[i] = ρL / Δt + ρL * uL[i+1] / Δx
            du_Lα[i] = 0.0
            dl_Lα[i] = -ρL * uL[i] / Δx
            b_Lα[i] = ρL * αL0[i] / Δt + Γ_L
            if uG[i] >= 0
                d_Gα[i] = ρG / Δt + ρG * uG[i+1] / Δx
                du_Gα[i] = 0.0
                dl_Gα[i] = -ρG * uG[i] / Δx
                b_Gα[i] = ρG * αG0[i] / Δt + Γ_G
            else
                d_Gα[i] = ρG / Δt - ρG * uG[i] / Δx
                du_Gα[i] = ρG * uG[i+1] / Δx
                dl_Gα[i] = 0.0
                b_Gα[i] = ρG * αG0[i] / Δt + Γ_G
            end
        end
        αLØ = TDMA(d_Lα, du_Lα, dl_Lα, b_Lα)
        αGØ = TDMA(d_Gα, du_Gα, dl_Gα, b_Gα)
        #
        # Correção das frações
        for j = 2:N
            αL[j] = αLØ[j] / (αLØ[j] + αGØ[j])
            αG[j] = 1.0 - αL[j]
        end
        #
        # Erros:  ##############################################################
        e_uL = norm(uL-uLØ)/norm(uLØ)
        e_uG = norm(uG-uGØ)/norm(uGØ)
        e_αL = norm(αL-αLØ)/norm(αLØ)
        e_αG = norm(αG-αGØ)/norm(αGØ)
        e_P = norm(P-PØ)/norm(PØ)
        #
        # Atualização dos Supostos -----------------------------------------
        uLØ[:] = uL[:]
        uGØ[:] = uG[:]
        αLØ[:] = αL[:]
        αGØ[:] = αG[:]
        PØ[:] = P[:]
    end
    #
    # Atualização Temporal -------------------------------------------------
    uL0[:] = uL[:]
    uG0[:] = uG[:]
    αL0[:] = αL[:]
    αG0[:] = αG[:]
    P0[:] = P[:]
end
#
df = DataFrame(Pr = P, UL = uL[1:N], UG = uG[1:N], AlfaL = αL, AlfaG = αG)
# CSV.write("./test.csv", df)
#plot 
fig = plot(0:Δx:Comp-Δx, αG)
savefig(fig, "./test.png")