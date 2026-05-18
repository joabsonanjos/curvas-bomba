analisar_bomba <- function(Q_novo, H_novo, N_novo, Hg_novo) {
  library(ggplot2)
  
  # 1. AJUSTES INICIAIS COM OS NOVOS VALORES ----
  curva <- data.frame(Q = Q_novo, H = H_novo, N = N_novo)
  Hg <- Hg_novo
  rotacoes <- c(1.0, 0.9, 0.8, 0.7)
  
  # 2. EXTRAÇÃO DE COEFICIENTES E MODELAGEM ----
  modelo_H <- lm(H ~ poly(Q, 2, raw = TRUE), data = curva)
  coefs_H  <- coef(modelo_H)
  c <- coefs_H[1]; b <- coefs_H[2]; a <- coefs_H[3]
  r2_H <- summary(modelo_H)$r.squared
  
  modelo_N <- lm(N ~ poly(Q, 2, raw = TRUE), data = curva)
  coefs_N  <- coef(modelo_N)
  f <- coefs_N[1]; e <- coefs_N[2]; d <- coefs_N[3]
  r2_N <- summary(modelo_N)$r.squared
  
  texto_H <- sprintf("Y = %.4fx² + %.4fx + %.4f\nR² = %.4f", a, b, c, r2_H)
  texto_N <- sprintf("Y = %.4fx² + %.4fx + %.4f\nR² = %.4f", d, e, f, r2_N)
  
  # Gráfico H vs Q Individual
  curva_QH <- ggplot(curva, aes(x = Q, y = H)) +
    geom_point(size = 2.5) +
    geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), color = "blue", se = FALSE) +
    annotate("text", x = max(curva$Q) * 0.6, y = max(curva$H) * 0.9, label = texto_H, hjust = 0, size = 4) +
    labs(title = "Ajuste Polinomial de Grau 2 (H vs Q)", x = "Q (m³/h)", y = "H (m)") +
    theme_minimal()
  
  # Gráfico N vs Q Individual
  curva_QN <- ggplot(curva, aes(x = Q, y = N)) +
    geom_point(size = 2.5) +
    geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), color = "red", se = FALSE) +
    annotate("text", x = max(curva$Q) * 0.1, y = max(curva$N) * 0.2, label = texto_N, hjust = 0, size = 4) +
    labs(title = "Ajuste Polinomial de Grau 2 (N vs Q)", x = "Q (m³/h)", y = "N (%)") +
    theme_minimal()
  
  # 3. DETERMINAÇÃO DO PONTO DE PROJETO E SISTEMA (K) ----
  Qreq <- -e / (2 * d)
  H1 <- a * Qreq^2 + b * Qreq + c
  N1 <- d * Qreq^2 + e * Qreq + f
  
  K <- (H1 - Hg) / Qreq^2
  
  curva$Q <- as.numeric(as.character(curva$Q))
  curva$HSIS <- Hg + (K * (curva$Q^2))
  
  porcentagens <- c(0.9, 0.8, 0.7)
  for(p in porcentagens) {
    nome_H <- paste0("H_", p*100, "pct")
    nome_N <- paste0("N_", p*100, "pct") 
    curva[[nome_H]] <- (a * (curva$Q^2)) + (b * p * curva$Q) + (c * (p^2))
    curva[[nome_N]] <- (d / (p^2) * (curva$Q^2)) + (e / p * curva$Q) + f
  }
  
  # 4. LEIS DE AFINIDADE E BHASKARA POR ROTAÇÃO ----
  tabela_completa <- data.frame(p = rotacoes)
  tabela_completa$a_p <- a
  tabela_completa$b_p <- b * tabela_completa$p
  tabela_completa$c_p <- c * (tabela_completa$p)^2
  
  calculo_bhaskara <- lapply(1:nrow(tabela_completa), function(i) {
    ap <- tabela_completa$a_p[i]; bp <- tabela_completa$b_p[i]; cp <- tabela_completa$c_p[i]
    A <- ap - K; B <- bp; C <- cp - Hg
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
  tabela_vazoes$H_op <- Hg + K * (tabela_vazoes$Q_op)^2
  tabela_vazoes$H_ph <- tabela_vazoes$H_op / (tabela_vazoes$p)^2
  
  # 6. CÁLCULO DO RENDIMENTO E COEFICIENTE K DA ROTAÇÃO ----
  tabela_vazoes$L_afin <- (d * (tabela_vazoes$Q_ph)^2) + (e * tabela_vazoes$Q_ph) + f
  tabela_vazoes$K_rot <- tabela_vazoes$H_ph / (tabela_vazoes$Q_ph)^2
  
  # 7. CRIAÇÃO DA TABELA: PARÁBOLAS HOMÓLOGAS (Q² x K_rot) ----
  tabela_rendimento <- data.frame(Q_ensaio = curva$Q)
  for(i in 1:nrow(tabela_vazoes)) {
    rotacao_decimal <- tabela_vazoes$p[i]
    K_atual <- tabela_vazoes$K_rot[i]
    nome_coluna <- sprintf("k_%.1f", rotacao_decimal)
    tabela_rendimento[[nome_coluna]] <- (curva$Q^2) * K_atual
  }
  
  # 8. GRAFICADOR COMBINADO DEFINITIVO ----
  grafico_final <- ggplot(curva, aes(x = Q)) +
    geom_smooth(aes(y = H, color = "H 100%"), method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, size = 0.8) +
    geom_smooth(aes(y = H_90pct, color = "H 90%"), method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, size = 0.8) +
    geom_smooth(aes(y = H_80pct, color = "H 80%"), method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, size = 0.8) +
    geom_smooth(aes(y = H_70pct, color = "H 70%"), method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, size = 0.8) +
    geom_smooth(aes(y = HSIS, color = "Sistema (HSIS)"), method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, size = 1) +
    geom_point(data = tabela_vazoes, aes(x = Q_op, y = H_op, color = "PF"), size = 3) +
    geom_smooth(aes(y = N, color = "Rendimento (Ensaio)"), method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, size = 1, linetype = "longdash") +
    geom_point(aes(y = N, color = "Rendimento (Ensaio)"), size = 2) + 
    geom_line(data = tabela_rendimento, aes(x = Q_ensaio, y = k_1.0, color = "K 100%"), linetype = "dashed", size = 0.8) +
    geom_line(data = tabela_rendimento, aes(x = Q_ensaio, y = k_0.9, color = "K 90%"), linetype = "dashed", size = 0.8) +
    geom_line(data = tabela_rendimento, aes(x = Q_ensaio, y = k_0.8, color = "K 80%"), linetype = "dashed", size = 0.8) +
    geom_line(data = tabela_rendimento, aes(x = Q_ensaio, y = k_0.7, color = "K 70%"), linetype = "dashed", size = 0.8) +
    scale_y_continuous(name = "Altura Manométrica / Parábolas Homólogas (m)", sec.axis = sec_axis(~ . * 1, name = "Rendimento (%)")) +
    scale_color_manual(values = c("H 100%"="blue","H 90%"="purple","H 80%"="green","H 70%"="orange","Sistema (HSIS)"="black","PF"="red","Rendimento (Ensaio)"="darkred","K 100%"="skyblue","K 90%"="plum","K 80%"="lightgreen","K 70%"="wheat")) +
    labs(title = "Curvas Características, Operação e Parábolas de Rendimento Constante", x = "Vazão Q (m³/h)", color = "Legenda") +
    theme_minimal()
  
  # 8.1 GRÁFICO ISOLADO RENDIMENTO ----
  curva_RR <- ggplot(curva, aes(x = Q)) +
    geom_smooth(aes(y = N, color = "N 100%"), method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, size = 1) +
    geom_smooth(aes(y = N_90pct, color = "N 90%"), method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, size = 1) +
    geom_smooth(aes(y = N_80pct, color = "N 80%"), method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, size = 1) +
    geom_smooth(aes(y = N_70pct, color = "N 70%"), method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, size = 1) +
    scale_color_manual(values = c("N 100%"="red","N 90%"="darkred","N 80%"="yellow","N 70%"="orange")) +
    labs(title = "Curvas de Rendimento por Patamar de Rotação", x = "Vazão Q (m³/h)", y = "Rendimento N (%)", color = "Patamares") +
    theme_minimal()
  
  write.csv2(tabela_vazoes, "tabela_vazoes_homologas.csv", row.names = FALSE)
  
  # Retorno estruturado em lista
  pacote_resultados <- list(
    tabela_ensaio_expandida = curva,
    tabela_pontos_operacao  = tabela_vazoes,
    tabela_parabolas_k      = tabela_rendimento,
    figura_ajuste_H         = curva_QH,
    figura_ajuste_N         = curva_QN,
    figura_combinada_final  = grafico_final,
    figura_rendimento_pure  = curva_RR
  )
  
  return(pacote_resultados)
}
