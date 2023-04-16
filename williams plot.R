
## El script recibe un dataset en CSV que contiene el CAS, el Tox y predictores. 
# Cargar librerías necesarias
library(caret)
library(pracma)

# Leer el archivo CSV (reemplaza 'tu_archivo.csv' por el nombre real de tu archivo)
data <- read.csv("ruta")

# Dividir el dataset en entrenamiento y prueba (en este caso son 32 para el training)
train_indices <- 1:32
test_indices <- 33:nrow(data)

# Extraer la columna Tox como la variable dependiente y- o LD50 o lo que sea
y_train <- data[train_indices, "Tox"]
y_test <- data[test_indices, "Tox"]

# Extraer las variables independientes (descriptores) en X
X <- data[, !(names(data) %in% c("CAS", "Tox"))]

# Dividir X en entrenamiento y prueba
X_train <- X[train_indices, ]
X_test <- X[test_indices, ]

# Ajusta un modelo de regresión lineal
train_data <- cbind(X_train, Tox = y_train)

# Ajustar un modelo de regresión lineal múltiple utilizando lm()
model <- lm(Tox ~ ., data = train_data)




## Función que calcula la diagonal de los hats. Debe recibir una matriz de predictores
#Sin nombres

hat_matrix = function(X){
  xtx = t(X) %*% X
  ixtx = solve(xtx) # Utiliza "solve()" en lugar de "inv()" para invertir la matriz
  return(X %*% ixtx %*% t(X))
}

## Esta es otra alternativa (Hace la misma mmda, pero usa una función nativa de R)
# hat_matrix = function(X){
#   xtx = t(X) %*% X
#   diag(xtx) = diag(xtx) + runif(n = length(diag(xtx)),min = 0.001,max = 0.002)
#   ixtx = inv(xtx)
#   return(X %*% (ixtx %*% t(X)))
# }

## Esta es una función que encontré de angy89/hyQSAR, pero tenía algunos errores.
## Ya se encuentra corregida 

williams_plot = function(X_train, X_test, y_train, y_test, model){
  beta <- coef(model)[-1] # Excluir el intercepto
  
  ## Aquí se calcula la matriz de apalancamientos
  
  H_train = hat_matrix(rbind(as.matrix(X_train), as.matrix(X_test)))
  
  ## Aquí se predicen los compuestos del training y del test. 

  y_pred_train = predict(model, X_train)
  y_pred_test = predict(model, X_test)
  
  #Aquí se calculan y se estandarizan los residuales.
  # Hay que tener en cuenta que el enfoque del estandarizamiento va más allá 
  ## Como lo veremos más adelante
  
  residual_train = abs(y_train - y_pred_train)
  residual_test = abs(y_test - y_pred_test)
  s_residual_train = (residual_train - mean(residual_train)) / sd(residual_train)
  s_residual_test = (residual_test - mean(residual_test)) / sd(residual_test)
  
  ## Aquí se dividen los apalancamientos anteriormente calculados en H_train
  
  n_train = nrow(X_train)
  n_test = nrow(X_test)
  leverage_train = H_train[1:n_train]
  leverage_test = H_train[(n_train + 1):(n_train + n_test)]
  
  # Calculamos los parámetros para la línea de advertencia
  
  p = length(beta) # features
  n = n_train # training compounds
  h_star = (3 * p) / n
  
  ## Aquí se utilizan los residuales estandarizados para operarlos con los leverages. 
  
  AD_train = 100 * (sum(leverage_train < h_star & abs(s_residual_train) < 3) / length(leverage_train))
  AD_test = 100 * (sum(leverage_test < h_star & abs(s_residual_test) < 3) / length(leverage_test))
  lev = c(leverage_train, leverage_test)
  res = c(s_residual_train, s_residual_test)
  col = c(rep("black", n_train), rep("red", n_test))
  
  data_to_plot = list(lev = lev, res = res, col = col, h_star = h_star)
  
  # Retorna los datos que necesitará la gráfica. 
  
  return(list(ADVal = c(AD_train, AD_test), DTP = data_to_plot))
}


# Aquí invocamos la función pasándole las variables que ya creamos anteriormente. 

data_to_plot <- williams_plot(X_train, X_test, y_train, y_test, model)

data_to_plot <- data_to_plot$DTP

## Función para crear el diagrama williams sin ggplot. 

plot_wp = function(data_to_plot){
  plot(x = data_to_plot$lev,y = data_to_plot$res,
       col =data_to_plot$col, ylim = c(min(-4.5,data_to_plot$res),max(3.5,data_to_plot$res)),
       xlim = c(0,max(max(data_to_plot$lev),data_to_plot$h_star)+0.1), xlab = "Leverage",ylab = "Standardized Residual")
  abline(v=data_to_plot$h_star)
  abline(h = 3, lty=3)
  abline(h = -3, lty=3)
  legend(x = "bottomright",legend = c("Train","Test"),fill = c("black","red"),bty = "n")
}
plot_wp( data_to_plot)



