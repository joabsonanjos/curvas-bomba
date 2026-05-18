
library(ggplot2)

#######################################################################
# 1. AJUSTES INICIAIS E CRIAÇÃO DO DATAFRAME ----
#######################################################################

Q <- c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45)
H <- c(76.5, 75.6, 74.3, 73, 70.4, 66.9, 60.8, 53, 42.1, 30)
N <- c(0, 22.2, 40.6, 55.1, 65.8, 72.6, 75.6, 74.8, 70.1, 61.5) 

Hg <- 15
rotacoes <- c(1.0, 0.9, 0.8, 0.7)

curva <- data.frame(Q = Q, H = H, N = N)
print("--- CURVA DO ENSAIO ---")
print(curva)

#######################################################################
# 1.1 GERAR OS GRÁFICOS
#######################################################################

# Alterado para %.4f em todos os coeficientes para exibir 4 casas decimais
texto_H <- sprintf("Y = %.4fx² + %.4fx + %.4f\nR² = %.4f", a, b, c, r2_H)
texto_N <- sprintf("Y = %.4fx² + %.4fx + %.4f\nR² = %.4f", d, e, f, r2_N)

# Gráfico H vs Q
curva_QH <- ggplot(curva, aes(x = Q, y = H)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), color = "blue", se = FALSE) +
  annotate("text", x = max(curva$Q) * 0.6, y = max(curva$H) * 0.9, label = texto_H, hjust = 0, size = 4) +
  labs(title = "Ajuste Polinomial de Grau 2 (H vs Q)", x = "Q (m³/h)", y = "H (m)") +
  theme_minimal()
print(curva_QH)

# Gráfico N vs Q
curva_QN <- ggplot(curva, aes(x = Q, y = N)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), color = "red", se = FALSE) +
  annotate("text", x = max(curva$Q) * 0.1, y = max(curva$N) * 0.2, label = texto_N, hjust = 0, size = 4) +
  labs(title = "Ajuste Polinomial de Grau 2 (N vs Q)", x = "Q (m³/h)", y = "N (%)") +
  theme_minimal()
print(curva_QN)

#######################################################################
# 2. EXTRAÇÃO DE COEFICIENTES E MODELAGEM ----
#######################################################################

# Coeficientes da Altura (H)
modelo_H <- lm(H ~ poly(Q, 2, raw = TRUE), data = curva)
coefs_H  <- coef(modelo_H)
c <- coefs_H[1] # Intercepto
b <- coefs_H[2] # Termo linear (X)
a <- coefs_H[3] # Termo quadrático (X²)
r2_H <- summary(modelo_H)$r.squared

# Coeficientes do Rendimento (N)
modelo_N <- lm(N ~ poly(Q, 2, raw = TRUE), data = curva)
coefs_N  <- coef(modelo_N)
f <- coefs_N[1] # Intercepto
e <- coefs_N[2] # Termo linear (X)
d <- coefs_N[3] # Termo quadrático (X²)
r2_N <- summary(modelo_N)$r.squared


#######################################################################
# 3. DETERMINAÇÃO DO PONTO DE PROJETO E SISTEMA (K) ----
#######################################################################

Qreq <- -e / (2 * d)
H1 <- a * Qreq^2 + b * Qreq + c
N1 <- d * Qreq^2 + e * Qreq + f

K <- (H1 - Hg) / Qreq^2

curva$Q <- as.numeric(as.character(curva$Q))
curva$HSIS <- Hg + (K * (curva$Q^2))

# Geração dos patamares de rotação para H e N na tabela inicial
porcentagens <- c(0.9, 0.8, 0.7)

for(p in porcentagens) {
  nome_H <- paste0("H_", p*100, "pct")
  nome_N <- paste0("N_", p*100, "pct") 
  
  # Afinidade aplicada à Altura Manométrica (H)
  curva[[nome_H]] <- (a * (curva$Q^2)) + (b * p * curva$Q) + (c * (p^2))
  
  # Afinidade aplicada ao Rendimento (N)
  curva[[nome_N]] <- (d / (p^2) * (curva$Q^2)) + (e / p * curva$Q) + f
}

print("--- TABELA INICIAL EXPANDIDA (HSIS, H E N ROTAÇÕES) ---")
print(curva)


#######################################################################
# 4. LEIS DE AFINIDADE E BHASKARA POR ROTAÇÃO ----
#######################################################################

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


#######################################################################
# 5. MONTAGEM DA TABELA FINAL DE VAZÕES E PRESSÕES ----
#######################################################################

colnames(tabela_completa) <- c("p", "a_p", "b_p", "c_p", "Q_op", "Q_ph")
tabela_vazoes <- tabela_completa[, c("p", "Q_op", "Q_ph")]

# Adicionando H_op e H_ph vetorizados
tabela_vazoes$H_op <- Hg + K * (tabela_vazoes$Q_op)^2
tabela_vazoes$H_ph <- tabela_vazoes$H_op / (tabela_vazoes$p)^2


#######################################################################
# 6. CÁLCULO DO RENDIMENTO E COEFICIENTE K DA ROTAÇÃO ----
#######################################################################

# Rendimento por afinidade usando a vazão homóloga (Q_ph)
tabela_vazoes$L_afin <- (d * (tabela_vazoes$Q_ph)^2) + (e * tabela_vazoes$Q_ph) + f

# Cálculo do K homólogo para cada patamar de rotação
tabela_vazoes$K_rot <- tabela_vazoes$H_ph / (tabela_vazoes$Q_ph)^2

print("--- TABELA DE OPERAÇÃO DINÂMICA (tabela_vazoes) ---")
print(tabela_vazoes)


#######################################################################
# 7. CRIAÇÃO DA TABELA: PARÁBOLAS HOMÓLOGAS (Q² x K_rot) ----
#######################################################################

tabela_rendimento <- data.frame(Q_ensaio = curva$Q)

for(i in 1:nrow(tabela_vazoes)) {
  rotacao_decimal <- tabela_vazoes$p[i]   # Pega o valor decimal (1.0, 0.9, 0.8...)
  K_atual <- tabela_vazoes$K_rot[i]       # Pega o K correspondente de cada patamar
  
  # Nome baseado no valor decimal: k_1.0, k_0.9, k_0.8...
  nome_coluna <- sprintf("k_%.1f", rotacao_decimal)
  
  # Vazão do ensaio ao quadrado MULTIPLICADA pelo K da rotação
  tabela_rendimento[[nome_coluna]] <- (curva$Q^2) * K_atual
}

print("--- TABELA: PARÁBOLAS HOMÓLOGAS (Q² * K) ---")
print(tabela_rendimento)


#######################################################################
# 8. GRAFICADOR COMBINADO DEFINITIVO ----
#######################################################################

grafico_final <- ggplot(curva, aes(x = Q)) +
  
  # A. Família de Curvas de Altura (H) - Linhas Sólidas (Eixo Esquerdo)
  geom_smooth(aes(y = H, color = "H 100%"), method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, size = 0.8) +
  geom_smooth(aes(y = H_90pct, color = "H 90%"), method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, size = 0.8) +
  geom_smooth(aes(y = H_80pct, color = "H 80%"), method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, size = 0.8) +
  geom_smooth(aes(y = H_70pct, color = "H 70%"), method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, size = 0.8) +
  
  # B. Curva do Sistema (HSIS) e Pontos de Operação reais (Bhaskara)
  geom_smooth(aes(y = HSIS, color = "Sistema (HSIS)"), method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, size = 1) +
  geom_point(data = tabela_vazoes, aes(x = Q_op, y = H_op, color = "PF"), size = 3) +
  
  # C. Curva de Rendimento do Ensaio (Eixo Direito - Linha Longa Tracejada)
  geom_smooth(aes(y = N, color = "Rendimento (Ensaio)"), method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, size = 1, linetype = "longdash") +
  geom_point(aes(y = N, color = "Rendimento (Ensaio)"), size = 2) + 
  
  # D. Parábolas de Rendimento Constante (Tracejadas Curtas baseadas na tabela_div_K)
  geom_line(data = tabela_div_K, aes(x = Q_ensaio, y = k_1.0, color = "K 100%"), linetype = "dashed", size = 0.8) +
  geom_line(data = tabela_div_K, aes(x = Q_ensaio, y = k_0.9, color = "K 90%"), linetype = "dashed", size = 0.8) +
  geom_line(data = tabela_div_K, aes(x = Q_ensaio, y = k_0.8, color = "K 80%"), linetype = "dashed", size = 0.8) +
  geom_line(data = tabela_div_K, aes(x = Q_ensaio, y = k_0.7, color = "K 70%"), linetype = "dashed", size = 0.8) +
  
  
  # E. Configuração dos Eixos (Esquerdo: Altura/Parábolas | Direito: Rendimento 1:1)
  scale_y_continuous(
    name = "Altura Manométrica / Parábolas Homólogas (m)",                          
    sec.axis = sec_axis(~ . * 1, name = "Rendimento (%)")     
  ) +
  
  # F. Mapeamento Manual de Cores da Legenda
  scale_color_manual(values = c(
    "H 100%"              = "blue",
    "H 90%"               = "purple",
    "H 80%"               = "green",
    "H 70%"               = "orange",
    "Sistema (HSIS)"      = "black",
    "PF"   = "red",
    "Rendimento (Ensaio)" = "darkred",
    # Cores das parábolas homólogas
    "K 100%"      = "skyblue",
    "K 90%"      = "plum",
    "K 80%"      = "lightgreen",
    "K 70%"      = "wheat"
  )) +
  
  labs(
    title = "Curvas Características, Operação e Parábolas de Rendimento Constante",
    x = "Vazão Q (m³/h)",
    color = "Legenda"
  ) +
  theme_minimal()

# Desenha o gráfico na aba Plots
print(grafico_final)

# =====================================================================
# 8.1 GRÁFICO ISOLADO: FAMÍLIA DE CURVAS DE RENDIMENTO (N vs Q)
# =====================================================================

curva_RR <- ggplot(curva, aes(x = Q)) +
  
  # 1. Plotagem de cada patamar de rotação do rendimento
  geom_smooth(aes(y = N, color = "N 100%"), method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, size = 1) +
  geom_smooth(aes(y = N_90pct, color = "N 90%"), method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, size = 1) +
  geom_smooth(aes(y = N_80pct, color = "N 80%"), method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, size = 1) +
  geom_smooth(aes(y = N_70pct, color = "N 70%"), method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, size = 1) +
  
  # Pontos originais do ensaio prático (opcional, apenas para o patamar de 100%)
  
  
  # 2. Configuração manual das cores da legenda
  scale_color_manual(values = c(
    "N 100%" = "red",
    "N 90%"  = "darkred",
    "N 80%"  = "yellow",
    "N 70%"  = "orange"
  )) +
  
  # 3. Títulos e rótulos dos eixos
  labs(
    title = "Curvas de Rendimento por Patamar de Rotação",
    x = "Vazão Q (m³/h)",
    y = "Rendimento N (%)",
    color = "Patamares"
  ) +
  theme_minimal()

plot(curva_RR)

#######################################################################
# 9. EXPORTAÇÃO DOS RESULTADOS ----
#######################################################################

write.csv2(tabela_vazoes, "tabela_vazoes_homologas.csv", row.names = FALSE)
