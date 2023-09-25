library(shiny)
library(ggplot2)
library(shinythemes)
library(rootSolve)
library(PlaneGeometry)

#Produces Ellipse of Insignificance Analysis (EOI) and Region of Attainable Redaction (ROAR) analysis 
#For details see: https://doi.org/10.7554/eLife.79573


#Enter a-d and alpha
a <- 700   
b <- 300
c <- 500
d <- 500
alpha <- 0.05

v <- qchisq(alpha, df=1, lower.tail=FALSE)
n <- a + b + c + d

chstat <- (n * ((a * d - b * c)^2)) / ((a + b) * (c + d) * (a + c) * (b + d)) 
pval <- pchisq(chstat, df=1, lower.tail=FALSE)

tstat_test <- chstat - v
vtest <- ifelse(tstat_test <= 0, 0, 1)

if (vtest == 0) {
  # output$output_text <- renderText({
  #  "Not significant before recoding!"
  output$table_output <- renderDataTable({
    table_dataerrmsg
  },options = list(searching = FALSE))
  output$plot1 <- renderPlot({})
  output$plot2 <- renderPlot({})
} else {
  # Calculations
  Ae <- (c + d) * ((c + d) * n + (a + b) * v)
  Be <- 2 * (a + b) * (c + d) * (n - v)
  Ce <- (a + b) * ((a + b) * n + (c + d) * v)
  De <- (c + d) * (2 * (b * c - a * d) * n + (a + b) * (b - a + d - c) * v)
  Ee <- (a + b) * (2 * (b * c - a * d) * n + (c + d) * (a - b + c - d) * v)
  Fe <- ((b * c - a * d)^2) * n - (a + b) * (a + c) * (b + d) * (c + d) * v
  
  # Find FECKUP vector points
  xp <- numeric()
  yp <- numeric()# Define the system of equations as functions
  # Define the equations
  
  # Define the equations
  # Define the equations
  equation1 <- function(xp, yp) {
    return ((2 * Ae * xp + Be * yp + De) * yp - xp * (Be * xp + 2 * Ce * yp + Ee))
  }
  
  equation2 <- function(xp, yp) {
    return (Ae * (xp^2) + Be * xp * yp + Ce * (yp^2) + De * xp + Ee * yp + Fe)
  }
  
  # Define a function that takes a vector of xp and yp and returns a vector of equations
  system_of_equations <- function(xy) {
    xp <- xy[1]
    yp <- xy[2]
    c(equation1(xp, yp), equation2(xp, yp))
  }
  
  # Number of solutions you want to find
  num_solutions <- 4
  
  # Initialize an empty list to store solutions
  solutions <- list()
  
  # Find all solutions
  attempts <- 0
  while (length(solutions) < num_solutions && attempts < 100) {
    # Generate random initial values
    initial_values <- c(runif(1, -n, n), runif(1, -n, n))
    
    # Use multiroot function to find roots
    result <- multiroot(system_of_equations, initial_values)
    
    # Calculate the distance between new solution and existing solutions
    distances <- sapply(solutions, function(sol) sum((sol - result$root)^2))
    
    # If the new solution is significantly different, store it
    if (all(distances > 1e-6)) {
      solutions <- c(solutions, list(result$root))
    }
    
    attempts <- attempts + 1
  }
  
  # Initialize a list to store the solution details
  solution_details <- list()
  
  # Display the solutions
  for (i in 1:length(solutions)) {
    #cat("Solution", i, ": xp =", solutions[[i]][1], ", yp =", solutions[[i]][2], "\n")
    sum_of_squares_i <- solutions[[i]][1]^2 + solutions[[i]][2]^2
    #cat("Sum of xp^2 + yp^2 for Solution", i, ":", sum_of_squares_i, "\n")
    
    # Store the solution details
    solution_details[[i]] <- list(xp = solutions[[i]][1], yp = solutions[[i]][2], sum_of_squares = sum_of_squares_i)
    
  }
  
  # Calculate the index with the smallest sum of squares
  dvmin <- which.min(sapply(solution_details, function(sol) sol$sum_of_squares))
  
  # Extract the values   
  xpv <- as.numeric(solution_details[[dvmin]][1])
  ypv <- as.numeric(solution_details[[dvmin]][2])
  FKUP2 <- as.numeric(solution_details[[dvmin]][3])
  FKUP <- sqrt(FKUP2)
  
  
  # Find xi point  
  A1 <- (c+d) * ((c+d) * n + (a+b) * v)
  B1 <- (c+d) * (2 * (b*c-a*d) * n + (a+b) * (b-a + d -c) * v)
  C1 <- ((b*c-a*d)^2) * n - (a+b) * (a+c) * (b+d) * (c+d) * v
  X1 <- (-B1 + sqrt(B1^2 - 4*A1*C1)) / (2*A1)
  X2 <- (-B1 - sqrt(B1^2 - 4*A1*C1)) / (2*A1)
  g <- c(X1, X2)
  gabs <- abs(g)
  xi <- g[which.min(gabs)]
  
  
  # Find yi point
  A2 <- (a + b) * ((a + b) * n + (c + d) * v)
  B2 <- (a + b) * (2 * (b * c - a * d) * n + (a - b + c - d) * (c + d) * v)
  C2 <- ((b * c - a * d)^2) * n - (a + b) * (a + c) * (b + d) * (c + d) * v
  Y1 <- (-B2 + sqrt(B2^2 - 4*A2*C2)) / (2*A2)
  Y2 <- (-B2 - sqrt(B2^2 - 4*A2*C2)) / (2*A2)
  g2 <- c(Y1, Y2)
  gabs2 <- abs(g)
  yi <- g2[which.min(gabs2)]
  
  #Older method of drawing ellipse
  xv <- seq(-1.05 * b, 1.05 * a, by = 0.01)
  By <- Be*xv + Ee;
  Cy <- Ae*(xv^2) + De*xv + Fe; 
  el1 <- (-By + sqrt((By^2) - 4 * Cy * Ce)) / (2 * Ce)
  el2 <- (-By - sqrt((By^2) - 4 * Cy * Ce)) / (2 * Ce)
  xplot <- c(min(xv), max(xv))
  yplot <- c(min(el2, na.rm = TRUE), max(el1, na.rm = TRUE))
  
  xtight <- sort(c(0-0.05*xi,1.05*xi))
  ytight <- sort(c(0-0.05*yi,1.05*yi))
  
  
  # Find all relevant terms
  dmin <- floor(abs(xpv) + abs(ypv))
  qvec <- c(1 - (a + b - abs(xi)) / (a + b), 
            1 - (c + d - abs(yi)) / (c + d), 
            1 - (n - abs(dmin)) / n)
  
  csize <- 1
  
  ratiocheck <- (a / (a + b)) / (c / (c + d))
  RRbound_upp <- exp(log(ratiocheck) + 1.96*sqrt((b/a)/(a+b) + (d/c)/(c+d)))
  RRbound_low <- exp(log(ratiocheck) - 1.96*sqrt((b/a)/(a+b) + (d/c)/(c+d)))
  
  if (ratiocheck > 1) {
    coeffx3y2 <- 1
    coeffx3y <- 2 * c
    coeffx3 <- c^2
    coeffx2y3 <- 1
    coeffx2y2 <- n + 2 * (b + c) - v
    coeffx2y <- 4 * b * c + c^2 - 2 * a * d + 2 * c * n - (a + 2 * c + d) * v
    coeffx2 <- c * (2 * b * c - 2 * a * d + c * n) - (a + c) * (c + d) * v
    coeffxy3 <- 2 * b
    coeffxy2 <- b^2 + 4 * b * c - 2 * a * d + 2 * b * n - (a + 2 * b + d) * v
    coeffxy <- 2 * (b + c) * (b * c - a * d) + 4 * b * c * n - 2 * a * d * n - (a + 2 * b + d) * (a + 2 * c + d) * v
    coeffx <- (b * c - a * d) * (b * c - a * d + 2 * c * n) - (a + c) * (a + 2 * b + d) * (c + d) * v
    coeffy3 <- b^2
    coeffy2 <- b * (-2 * a * d + b * (2 * c + n)) - (a + b) * (b + d) * v
    coeffy <- (b * c - a * d) * (-a * d + b * (c + 2 * n)) - (a + b) * (b + d) * (a + 2 * c + d) * v
    coeffc <- (b * c - a * d)^2 * n - (a + b) * (a + c) * (b + d) * (c + d) * v
    
    gfull <- function(x, y) {
      return(coeffx3y2 * (x^3 * y^2) + coeffx3y * (x^3 * y) + coeffx3 * (x^3) +
               coeffx2y3 * (x^2 * y^3) + coeffx2y2 * (x^2 * y^2) + coeffx2y * (x^2 * y) +
               coeffx2 * (x^2) + coeffxy3 * (x * y^3) + coeffxy2 * (x * y^2) +
               coeffxy * (x * y) + coeffx * (x) + coeffy3 * (y^3) + coeffy2 * (y^2) +
               coeffy * (y) + coeffc)
    }
    
    disfun <- function(x, y) {
      return(sqrt(x^2 + y^2))
    }
    
    
    dgdx <- function(x, y) {
      return(3*coeffx3y2 * (x^2 * y^2) + 3*coeffx3y * (x^2 * y) + 3*coeffx3 * (x^2) +
               2*coeffx2y3 * (x * y^3) + 2*coeffx2y2 * (x * y^2) + 2*coeffx2y * (x * y) +
               2*coeffx2 * (x) + coeffxy3 * (y^3) + coeffxy2 * ( y^2) +
               coeffxy * (y) + coeffx * (1)) 
    }
    
    dgdy <-  function(x, y) {
      return(2*coeffx3y2 * (x^3 * y) + coeffx3y * (x^3 * 1) + 0 +
               3*coeffx2y3 * (x^2 * y^2) + 2*coeffx2y2 * (x^2 * y) + coeffx2y * (x^2 * 1) +
               0 + 3*coeffxy3 * (x * y^2) + 2*coeffxy2 * (x * y) +
               coeffxy * (x * 1) + 0 + 3*coeffy3 * (y^2) + 2*coeffy2 * (y) +
               coeffy * (1) + 0)
    }
    
    
    
    
    
    dDdx <- function(x, y) {
      return(x/sqrt(x^2 + y^2))
    }
    
    
    dDdy <- function(x, y) {
      return(y/sqrt(x^2 + y^2))
    }
    
    eqsimp <- function(x, y) {
      return(y * dgdx(x, y) - x * dgdy(x, y))
    }
    
    
    # Define a function that takes a vector of xp and yp and returns a vector of equations
    system_of_equations <- function(xy) {
      x <- xy[1]
      y <- xy[2]
      c(gfull(x, y), eqsimp(x, y))
    }
    
    # Number of solutions you want to find
    num_solutions <- 25
    
    # Initialize an empty list to store solutions
    solutions <- list()
    
    # Find all solutions
    attempts <- 0
    while (length(solutions) < num_solutions && attempts < 500) {
      # Generate random initial values
      initial_values <- c(runif(1, -0.1*n, n), runif(1, -0.1*n, n))
      
      # Use multiroot function to find roots
      result <- multiroot(system_of_equations, initial_values)
      
      # Calculate if root is + 
      rsign1 <- result$root[[1]]
      rsign2 <- result$root[[2]]
      
      # Calculate the distance between new solution and existing solutions
      distances <- sapply(solutions, function(sol) sum((sol - result$root)^2))
      
      # If the new solution is significantly different, store it
      if (all(distances > 1e-6) && rsign1 > 0 && rsign2 > 0) {
        solutions <- c(solutions, list(result$root))
      }
      
      attempts <- attempts + 1
    }
    
    # Initialize a list to store the solution details
    solution_details <- list()
    
    # Display the solutions
    for (i in 1:length(solutions)) {
      #cat("Solution", i, ": x =", solutions[[i]][1], ", y =", solutions[[i]][2], "\n")
      sum_of_squares_i <- solutions[[i]][1]^2 + solutions[[i]][2]^2
      #cat("Sum of x^2 + y^2 for Solution", i, ":", sum_of_squares_i, "\n")
      
      # Store the solution details
      solution_details[[i]] <- list(xp = solutions[[i]][1], yp = solutions[[i]][2], sum_of_squares = sum_of_squares_i)
      
    }
    
    # Calculate the index with the smallest sum of squares
    dvmin <- which.min(sapply(solution_details, function(sol) sol$sum_of_squares))
    
    # Extract the values   
    xe <- as.numeric(solution_details[[dvmin]][1])
    ye <- as.numeric(solution_details[[dvmin]][2])
    FOCK <- sqrt(as.numeric(solution_details[[dvmin]][3]))
    
    #now we solve for xc and yc. first define...
    
    
    gyfull <- function(y) {
      return(0 + 0 + 0 +
               0 + 0 + 0 +
               0 + 0 + 0 +
               0 + 0 + coeffy3 * (y^3) + coeffy2 * (y^2) +
               coeffy * (y) + coeffc)
    }
    
    gxfull <- function(x) {
      return(0 + 0 + coeffx3 * (x^3) +
               0 + 0 + 0 +
               coeffx2 * (x^2) + 0 + 0 +
               0 + coeffx * (x) + 0 + 0 +
               0 + coeffc)
    }
    
    
    #search for 2 roots around xe,ye
    # Number of solutions you want to find
    num_solutions <- 3
    
    # Initialize an empty list to store solutions
    solutionsxc <- list()
    solutionsyc <- list()
    
    # Find all solutions
    attempts <- 0
    while (length(solutions) < num_solutions && attempts < 100) {
      
      # Use multiroot function to find roots
      resultxc <- multiroot(gxfull, xe)
      resultyc <- multiroot(gyfull, ye)
      
      # Calculate if root is + 
      rsignx <- resultxc$root
      rsigny <- resultyc$root
      
      # Calculate the distance between new solution and existing solutions
      distancesxc <- sapply(solutionsxc, function(sol) sum((sol - resultxc$root)^2))
      distancesyc <- sapply(solutionsyc, function(sol) sum((sol - resultyc$root)^2))
      
      # If the new solution is significantly different, store it
      if (all(distancesxc > 1e-6) && rsignx > 0) {
        solutionsxc <- c(solutionsxc, list(resultxc$root))
      }
      
      # If the new solution is significantly different, store it
      if (all(distancesyc > 1e-6) && rsigny > 0) {
        solutionsyc <- c(solutionsyc, list(resultyc$root))
      }
      
      
      attempts <- attempts + 1
    }
    
    yc <- min(as.numeric(solutionsyc))
    xc <- min(as.numeric(solutionsxc))
    
    
    
  } else {
    coeffx3y2 <- 1
    coeffx3y <- 2*d
    coeffx3 <- d^2
    coeffx2y3 <- 1
    coeffx2y2 <- n + 2*(a + d) - v
    coeffx2y <- -2*b*c + d*(4*a + d + 2*n) - v*(b + c + 2*d) 
    coeffx2 <- d*(-2 *b* c + d*(2 *a + n)) - (b + d)*(c + d)*v
    coeffxy3 <- 2*a
    coeffxy2 <- -2*b*c + a*(a + 4*d + 2*n) - (2*a + b + c)*v
    coeffxy <- 2*(a + d)*(-b*c + a*d) - 2*b*c*n + 4*a*d*n - (2*a + b + c)*(b + c + 2*d)*v
    coeffx <- (b*c - a*d)*(b*c - d*(a + 2*n)) - (2*a + b + c)*(b + d)*(c + d)*v
    coeffy3 <- a^2
    coeffy2 <- a*(-2*b*c + a*(2*d + n)) - (a + b)*(a + c)*v
    coeffy <- (b*c - a*d)*(b*c - a*(d + 2*n)) - (a + b)*(a + c)*(b + c + 2*d)*v
    coeffc <- ((b*c - a*d)^2)*n - (a + b)*(a + c)*(b + d)*(c + d)*v
    
    gfull <- function(x, y) {
      return(coeffx3y2 * (x^3 * y^2) + coeffx3y * (x^3 * y) + coeffx3 * (x^3) +
               coeffx2y3 * (x^2 * y^3) + coeffx2y2 * (x^2 * y^2) + coeffx2y * (x^2 * y) +
               coeffx2 * (x^2) + coeffxy3 * (x * y^3) + coeffxy2 * (x * y^2) +
               coeffxy * (x * y) + coeffx * (x) + coeffy3 * (y^3) + coeffy2 * (y^2) +
               coeffy * (y) + coeffc)
    }
    
    disfun <- function(x, y) {
      return(sqrt(x^2 + y^2))
    }
    
    
    dgdx <- function(x, y) {
      return(3*coeffx3y2 * (x^2 * y^2) + 3*coeffx3y * (x^2 * y) + 3*coeffx3 * (x^2) +
               2*coeffx2y3 * (x * y^3) + 2*coeffx2y2 * (x * y^2) + 2*coeffx2y * (x * y) +
               2*coeffx2 * (x) + coeffxy3 * (y^3) + coeffxy2 * ( y^2) +
               coeffxy * (y) + coeffx * (1)) 
    }
    
    dgdy <-  function(x, y) {
      return(2*coeffx3y2 * (x^3 * y) + coeffx3y * (x^3 * 1) + 0 +
               3*coeffx2y3 * (x^2 * y^2) + 2*coeffx2y2 * (x^2 * y) + coeffx2y * (x^2 * 1) +
               0 + 3*coeffxy3 * (x * y^2) + 2*coeffxy2 * (x * y) +
               coeffxy * (x * 1) + 0 + 3*coeffy3 * (y^2) + 2*coeffy2 * (y) +
               coeffy * (1) + 0)
    }
    
    
    
    
    
    dDdx <- function(x, y) {
      return(x/sqrt(x^2 + y^2))
    }
    
    
    dDdy <- function(x, y) {
      return(y/sqrt(x^2 + y^2))
    }
    
    eqsimp <- function(x, y) {
      return(y * dgdx(x, y) - x * dgdy(x, y))
    }
    
    
    # Define a function that takes a vector of xp and yp and returns a vector of equations
    system_of_equations <- function(xy) {
      x <- xy[1]
      y <- xy[2]
      c(gfull(x, y), eqsimp(x, y))
    }
    
    # Number of solutions you want to find
    num_solutions <- 25
    
    # Initialize an empty list to store solutions
    solutions <- list()
    
    # Find all solutions
    attempts <- 0
    while (length(solutions) < num_solutions && attempts < 500) {
      # Generate random initial values
      initial_values <- c(runif(1, -0.1*n, n), runif(1, -0.1*n, n))
      
      # Use multiroot function to find roots
      result <- multiroot(system_of_equations, initial_values)
      
      # Calculate if root is + 
      rsign1 <- result$root[[1]]
      rsign2 <- result$root[[2]]
      
      # Calculate the distance between new solution and existing solutions
      distances <- sapply(solutions, function(sol) sum((sol - result$root)^2))
      
      # If the new solution is significantly different, store it
      if (all(distances > 1e-6) && rsign1 > 0 && rsign2 > 0) {
        solutions <- c(solutions, list(result$root))
      }
      
      attempts <- attempts + 1
    }
    
    # Initialize a list to store the solution details
    solution_details <- list()
    
    # Display the solutions
    for (i in 1:length(solutions)) {
      #cat("Solution", i, ": x =", solutions[[i]][1], ", y =", solutions[[i]][2], "\n")
      sum_of_squares_i <- solutions[[i]][1]^2 + solutions[[i]][2]^2
      #cat("Sum of x^2 + y^2 for Solution", i, ":", sum_of_squares_i, "\n")
      
      # Store the solution details
      solution_details[[i]] <- list(xp = solutions[[i]][1], yp = solutions[[i]][2], sum_of_squares = sum_of_squares_i)
      
    }
    
    # Calculate the index with the smallest sum of squares
    dvmin <- which.min(sapply(solution_details, function(sol) sol$sum_of_squares))
    
    # Extract the values   
    xe <- as.numeric(solution_details[[dvmin]][1])
    ye <- as.numeric(solution_details[[dvmin]][2])
    FOCK <- sqrt(as.numeric(solution_details[[dvmin]][3]))
    
    
    
    
    #now we solve for xc and yc. first define...
    
    
    gyfull <- function(y) {
      return(0 + 0 + 0 +
               0 + 0 + 0 +
               0 + 0 + 0 +
               0 + 0 + coeffy3 * (y^3) + coeffy2 * (y^2) +
               coeffy * (y) + coeffc)
    }
    
    gxfull <- function(x) {
      return(0 + 0 + coeffx3 * (x^3) +
               0 + 0 + 0 +
               coeffx2 * (x^2) + 0 + 0 +
               0 + coeffx * (x) + 0 + 0 +
               0 + coeffc)
    }
    
    
    #search for 2 roots around xe,ye
    # Number of solutions you want to find
    num_solutions <- 3
    
    # Initialize an empty list to store solutions
    solutionsxc <- list()
    solutionsyc <- list()
    
    # Find all solutions
    attempts <- 0
    while (length(solutions) < num_solutions && attempts < 100) {
      
      # Use multiroot function to find roots
      resultxc <- multiroot(gxfull, xe)
      resultyc <- multiroot(gyfull, ye)
      
      # Calculate if root is + 
      rsignx <- resultxc$root
      rsigny <- resultyc$root
      
      # Calculate the distance between new solution and existing solutions
      distancesxc <- sapply(solutionsxc, function(sol) sum((sol - resultxc$root)^2))
      distancesyc <- sapply(solutionsyc, function(sol) sum((sol - resultyc$root)^2))
      
      # If the new solution is significantly different, store it
      if (all(distancesxc > 1e-6) && rsignx > 0) {
        solutionsxc <- c(solutionsxc, list(resultxc$root))
      }
      
      # If the new solution is significantly different, store it
      if (all(distancesyc > 1e-6) && rsigny > 0) {
        solutionsyc <- c(solutionsyc, list(resultyc$root))
      }
      
      
      attempts <- attempts + 1
    }
    
    yc <- min(as.numeric(solutionsyc))
    xc <- min(as.numeric(solutionsxc))
    
  }
  
  #tempplotfunc
  # Create a grid of x and y values
  xgrid <- seq(0, 1.05*xc, by = 1.05*xc/250)
  ygrid <- seq(0, 1.05*yc, by = 1.05*yc/250)
  fullgrid <- expand.grid(x = xgrid, y = ygrid)
  
  z <- matrix(0, nrow = length(xgrid), ncol = length(ygrid))
  for (i in 1:length(xgrid)) {
    for (j in 1:length(ygrid)) {
      if (gfull(xgrid[i], ygrid[j]) > 0) {
        z[i, j] <- 0
      } else {
        z[i, j] <- 1
      }
    }
  }
  
  #stats on FOCK vector
  rmin <- floor(xe + ye)
  rhoe <- 1 - (a+b)/(a+b+xc)
  rhoc <- 1 - (c+d)/(c+d+yc)
  rhoall <- 1 - (n)/(n + rmin)
  resulthrr <- paste("Reported Risk Ratio (Experimental - Control) :", ratiocheck, " (",  RRbound_low, " - ", RRbound_upp, ")", sep = "")
  
  nrows <- length(xgrid)
  ncols <- length(ygrid)
  
  angche <- 180*(atan(ye / xe))/pi
  
  
  # Create a data frame for ggplot
  df <- data.frame(
    x = rep(xgrid, each = ncols),
    y = rep(ygrid, times = nrows),
    z = as.vector(z)
  )
  
  
  ell <- EllipseFromEquation(A = Ae, B = Be, C = Ce, D = De, E = Ee, F = Fe)
  origincirc <- Circle$new(c(0,0), csize)
  xicirc <- Circle$new(c(xi,0), csize)
  yicirc <- Circle$new(c(0,yi), csize)
  fkcirc <- Circle$new(c(xpv,ypv), csize)
  box <- ell$boundingbox()
  plot(NULL, xlim = xtight, ylim = ytight, xlab = "Experimental group re-coding (x)", ylab = "Control group re-coding (y)", main = "Ellipse of insignificance (region of interest)")
  draw(ell, col = "lightblue", border = "purple", lwd = 3)
  draw(origincirc, col = rgb(0, 0, 0, alpha = 0.5), border = "black", lwd = 1)
  draw(xicirc, col = rgb(0, 0, 0, alpha = 0.5), border = "black", lwd = 1)
  #draw(yicirc, col = "white", border = "black", lwd = 1)
  draw(yicirc, col = rgb(0, 0, 0, alpha = 0.5), border = "black", lwd = 1)
  draw(fkcirc, col = rgb(0, 0, 0, alpha = 0.5), border = "black", lwd = 1)
  segments(x0 = 0, y0 = 0, x1 = xi, y1 = 0, col = "darkgreen", lty = "dashed") 
  segments(x0 = 0, y0 = 0, x1 = 0, y1 = yi, col = "darkgreen", lty = "longdash") 
  segments(x0 = 0, y0 = 0, x1 = xpv, y1 = ypv, col = "red", lty = "dotted", lwd = 3) 
  text(c(xpv/2),c(ypv/2),labels=c("FECKUP"))
  text(c(0),c(yi/2),labels=c("|yi|"))
  text(c(xi/2),c(0),labels=c("|xi|"))
  
 ggplot(df, aes(x = x, y = y)) +
              geom_tile(aes(fill = factor(z))) + guides(fill = FALSE) +
              scale_fill_manual(values = c("white", "lightblue")) +
              labs(x = "Experimental arm redaction", y = "Control arm redaction", fill = "Z Values", title="Region of Attainable redaction (ROAR)") +
              theme_classic() + theme(plot.title = element_text(color = "black", size = 15, face = "bold")) +
              annotate("segment", x = 0, y = 0, xend = 0, yend = yc, colour = "darkgreen", linetype = "dashed", size = 2) +
              annotate("segment", x = 0, y = 0, xend = xc, yend = 0, colour = "blue", linetype = "dashed", size = 2) +
              annotate("segment", x = 0, y = 0, xend = xe, yend = ye, colour = "red", size = 2) +  
              annotate("segment", x = 0, y = 0, xend = xe, yend = 0, colour = "orange", linetype = "dotted", size = 1) +
              annotate("segment", x =xe, y = 0, xend = xe, yend = ye, colour = "orange", linetype = "dotted", size = 1) +
              annotate("text",x = xe / 2, y = 1, label = "|xe|", size = 3, color = "black") +
              annotate("text",x = xe + 1, y = ye/2, label = "|ye|", size = 3, color = "black") +
              annotate("text",x = xc/2, y = 1, label = "|xc|", size = 3, color = "black") +
              annotate("text",x = 1, y = yc/2, label = "|yc|", size = 3, color = "black") +
              annotate("text",x = xe / 2, y = ye / 2, label = "FOCK", size = 4, color = "black", angle = angche)
  
  
  }
    