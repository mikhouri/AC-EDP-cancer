# Arquivo das funções da simulação em automato.jl

#Calcula o coeficinte que será usado na aproximação numérica
function calcula_coef(C, alfa, lambda, h, n)
    v = zeros(n+1,n+1)
    sign = 0 # célula normal -> função da posição
    sigc = 0 # célula cancerígena -> função da posição
    for i = 1:n+1
      for j = 2:n
        # célula normal  
        if C[i,j] == 0
          sigc = 0
          sign = 1
          v[i,j] = (h^2)*(alfa^2)*(sign + lambda*sigc)
        # célula cancerígena
        elseif C[i,j] > 0
          sigc = C[i,j]
          sign = 1
          v[i,j] = (h^2)*(alfa^2)*(sign + lambda*sigc)
        # célula morta (não tem consumo)
        else
          v[i,j] = 0
        end
      end
    end
    return v
end

# Calcula a solução analítica das EDP's
function analitica(n, alfa, sign, sigc, lambda)
    # Coeficientes
    c = (alfa^2)*(sign + lambda*sigc)
    A = exp(-sqrt(c)*n)/(exp(sqrt(c)*n)+exp(-sqrt(c)*n))
    B = 1-A
    # Solução geral
    phi = zeros(n+1)
    x = linspace(0,n+1,n+1)
    phi = A*exp.(sqrt(c)*x) + B*exp.(-sqrt(c)*x)
    return phi
end

# Calcula a aproximação numérica através do SOR
function sor(w, v, n, alfa, sign, sigc, lambda)
    U = zeros(n+1,n+1)
    phi = analitica(n, alfa, sign, sigc, lambda)    
    for i=1:n+1 U[i,:] = phi'  end
    # Fronteira inicial
    U[:,1] = 1
    # r: número de iterações
    r = 250
    for k = 1:r
        for i = 1:n+1
          # Fronteira superior
          if i == 1
            for j = 2:n
               U[i,j] = (1-w)*U[i,j] + w*(1/(4+v[i,j]))*(U[n+1,j]+U[i+1,j]+U[i,j-1]+U[i,j+1])
            end
          # Fronteira inferior  
          elseif i == n+1
            for j = 2:n
               U[i,j] = (1-w)*U[i,j] + w*(1/(4+v[i,j]))*(U[i-1,j]+U[1,j]+U[i,j-1]+U[i,j+1])
            end   
          # Restante da Matriz  
          else
            for j = 2:n
               U[i,j] = (1-w)*U[i,j] + w*(1/(4+v[i,j]))*(U[i-1,j]+U[i+1,j]+U[i,j-1]+U[i,j+1])
            end
          end
        end      
    end
    return U
end

# Encontra o número total de células cancerígenas no grid do automato
function conta_cancer(C,n)::Int64
    k = 0
    for i = 2:n
        for j = 2:n
            if C[i,j] > 0
                k = k + C[i,j]
            end
        end
    end
    return k 
end    

# Cria vetores para armazenar os índices das células cancerígenas
function ordena(C, q, n, b)::Array{Int64,1}
  I = Array{Int64,1}(q)
  J = Array{Int64,1}(q)
  o = 1
  for i = 2:n
    for j= 2:n
      if C[i,j] == 1       
        I[o] = i
        J[o] = j
        o = o+1
      elseif C[i,j] > 1
        c = Int(C[i,j])
        I[o:o+c-1] = i
        J[o:o+c-1] = j
        o = o+c
      end
    end
  end
  if b == 1   
    return I
  elseif b == 2  
    return J
  end
end

# Função que avalia a vizinhança da célula que foi marcada para se mover e executa a migração caso ela seja factível
function move(C,i,j)
    viz = [i+1 i-1 j+1 j-1]
    idx = rand(1:4)
    # Célula vizinha NÃO é cancerígena
    if idx <=2 && C[viz[idx], j] < 1
        C[viz[idx],j] = 1
        C[i,j] = C[i,j] - 1
    elseif idx > 2 && C[i, viz[idx]] < 1
        C[i, viz[idx]] = 1
        C[i,j] = C[i,j] - 1
    end
    # Célula vizinha É cancerígena: não se move 
    return C
end
