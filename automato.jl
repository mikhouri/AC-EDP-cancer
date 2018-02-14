# Rotina principal, consiste na evolução do autômato celular

include("/.../funcoes.jl")

n = 500
h = 1
lambdaN = 50
lambdaM = 100
alfa = 2/n

C = zeros(n+1,n+1) # Grid do AC

v = calcula_coef(C, alfa, lambdaN, h, n)
f = calcula_coef(C, alfa, lambdaM, h, n)

N = sor(1.7, v, n, alfa, 1, 0, lambdaN)
M = sor(1.7, f, n, alfa, 1, 0, lambdaM)

# Parâmetros de forma das probabillidades
theta_div = 0.3
theta_mov = 2
theta_del = 0.045

C[250,300] = 1

# 1 == mov
# 2 == div
# 3 == del

q = conta_cancer(C,n)
G = q

for t = 1:750
    
    I = ordena(C,q,n,1)
    J = ordena(C,q,n,2)
    
    while q != 0
	 
        # Sorteia célula cancerígena
        a = rand(1:q)::Int64
        i = I[a]
	      j = J[a]
	
        # Sorteia ação
        # quando theta_mov = infty sorteamos valores do vetor [2 3]
        act = rand(1:3) 
	
        if act == 1
            p_mov = 1-exp(-C[i,j]*(N[i,j]/theta_mov)^2)
            r = rand()
            if r <= p_mov
                C = move(C,i,j)
            end
	
        elseif act == 2           
            p_div = 1-exp(-(N[i,j]/(C[i,j]*theta_div))^2)     
            r = rand()  
            if r <= p_div
                C[i,j] = C[i,j] + 1   
		            G = G + 1
                C = move(C,i,j)     
            end
		
        elseif act == 3
            p_del = exp(-(M[i,j]/(C[i,j]*theta_del))^2)
            r = rand()
            if r <= p_del		
		            if C[i,j] == 1
            		    C[i,j] = -1
		            elseif C[i,j] > 1
			              C[i,j] = C[i,j] -1
		            end
		            G = G - 1
            end
        end    
	
	      # Atualiza parâmetros para próxima célula: remove do vetor de índices a células que já foi selecionada
      	I[a] = -10
        J[a] = -10
        o = 1
        IT = Array{Int64,1}(q-1)
        JT = Array{Int64,1}(q-1)
        for m = 1:q

            if I[m] > 0 && J[m] > 0
                IT[o] = I[m]
                JT[o] = J[m]
                o = o+1
            end
        end

        I = Array{Int64,1}(q-1)
        J = Array{Int64,1}(q-1)
        I = IT
        J = JT
        q = q - 1
    end

    # Atualiza parâmetros para próxima iteração
    v = calcula_coef(C, alfa, lambdaN, h, n)
    f = calcula_coef(C, alfa, lambdaM, h, n)

    N = sor(1.7, v, n, alfa, 1, 0, lambdaN)
    M = sor(1.7, f, n, alfa, 1, 0, lambdaM)
    q = G
end

# Gera a imagem final do tumor
using Plots
P = figura(C)
display(P)
