simular_bomba <- function(Q, H, N, Hg = 3) {

library(ggplot2)
library(ggpmisc)

###################
# 1. AJUSTES INICIAIS E GRÁFICOS ----

# Criando a tabela (Data Frame) corretamente com colunas nomeadas
curva <- data.frame(Q = Q, H = H, N = N)
print(curva)

###################

# Gráfico H vs Q com Equação
df_plot <- data.frame(Q = Q, H = H)
ggplot(df_plot, aes(x = Q, y = H)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), color = "blue") +
  stat_poly_eq(formula = y ~ poly(x, 2, raw = TRUE), 
               aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "*\", \"*")), 
               parse = TRUE, label.x = "right", label.y = "top") +
  labs(title = "Ajuste Polinomial de Grau 2 (H vs Q)", x = "Vazão (Q)", y = "Altura (H)") +
  theme_minimal()

# Gráfico N vs Q com Equação
ggplot(curva, aes(x = Q, y = N)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), color = "red") +
  stat_poly_eq(formula = y ~ poly(x, 2, raw = TRUE), 
               aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "*\", \"*")), 
               parse = TRUE) +
  labs(title = "Ajuste Polinomial de Grau 2 (N vs Q)", x = "Vazão (Q)", y = "Rendimento (N)") +
  theme_minimal()


# 2. EXTRAÇÃO DE COEFICIENTES E MODELAGEM DELA ----
# Coeficientes da Altura (H)
modelo_H <- lm(H ~ poly(Q, 2, raw = TRUE), data = curva)
summary(modelo_H) # Corrigido aqui
coefs_H <- coef(modelo_H)
c <- coefs_H[1]
b <- coefs_H[2]
a <- coefs_H[3]

# Coeficientes do Rendimento (N)
modelo_N <- lm(N ~ poly(Q, 2, raw = TRUE), data = curva)
summary(modelo_N)
coefs_N <- coef(modelo_N)
f <- coefs_N[1]
e <- coefs_N[2]
d <- coefs_N[3]


# 3. DETERMINAÇÃO DO PONTO DE PROJETO E SISTEMA (K) ----
Qreq <- -e / (2 * d)
H1 <- a * Qreq^2 + b * Qreq + c
N1 <- d * Qreq^2 + e * Qreq + f


K <- (H1 - Hg) / Qreq^2

curva$Q <- as.numeric(as.character(curva$Q))
curva$HSIS <- Hg + (K * (curva$Q^2))

print(curva$HSIS)

# 4. LEIS DEHSIS# 4. LEIS DE AFINIDADE E BHASKARA POR ROTAÇÃO ----

tabela_completa <- data.frame(p = rotacoes)

tabela_completa$a_p <- a
tabela_completa$b_p <- b * tabela_completa$p
tabela_completa$c_p <- c * (tabela_completa$p)^2

# Execução do loop de Bhaskara para achar cruzamentos reais
calculo_bhaskara <- lapply(1:nrow(tabela_completa), function(i) {
  ap <- tabela_completa$a_p[i]
  bp <- tabela_completa$b_p[i]
  cp <- tabela_completa$c_p[i]
  
  A <- ap - K
  B <- bp
  C <- cp - Hg
  
  delta <- (B^2) - (4 * A * C)
  q_op  <- (-B - sqrt(delta)) / (2 * A)
  q_ph  <- q_op / tabela_completa$p[i]
  
  return(c(Q_op = q_op, Q_ph = q_ph))
})

resultados_finais <- do.call(rbind, calculo_bhaskara)
tabela_completa <- cbind(tabela_completa, resultados_finais)


# 5. MONTAGEM DA TABELA FINAL DE VAZÕES E PRESSÕES ----
colnames(tabela_completa) <- c("p", "a_p", "b_p", "c_p", "Q_op", "Q_ph")
tabela_vazoes <- tabela_completa[, c("p", "Q_op", "Q_ph")]

# Adicionando H_op e H_ph vetorizados (sem o [i] fantasma)
tabela_vazoes$H_op <- Hg + K * (tabela_vazoes$Q_op)^2
tabela_vazoes$H_ph <- tabela_vazoes$H_op / (tabela_vazoes$p)^2


# 6. CÁLCULO DO RENDIMENTO POR AFINIDADE (FALTA DO SCRIPT ORIGINAL) ----
# Usando a vazão homóloga (Q_ph) na curva de rendimento de 100% (d, e, f)
tabela_vazoes$L_afin <- (d * (tabela_vazoes$Q_ph)^2) + (e * tabela_vazoes$Q_ph) + f


# 7. EXPORTAÇÃO DOS RESULTADOS ----

return(tabela_vazoes)
}
