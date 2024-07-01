calc.J.res <- function(w, I, J) { # I <- c(1,4)     # J <- c(2,3,5)
  J.res <- matrix(nrow = nrow(w), ncol = length(J))
  if (length(I) == 0) { # Check for root node
    for (j in J) {
      J.res[,which(J == j)] <- lm(w[,j]~1, data = w[I])$residuals
    }
  } else if (length(J) == 0) { # Check for leaf node
    return(NULL)
  } else { # Intermediate node
    for (j in J) {
      J.res[,which(J == j)] <- lm(w[,j]~., data = w[I])$residuals
    }
  }
  return(data.frame(J.res))
}

f=function(p,sr,lx2,ux2,ly2,uy2,lxy,uxy) {
  tol=.00001
  lx=sqrt(lx2)
  ux=sqrt(ux2)
  ly=sqrt(ly2)
  uy=sqrt(uy2)
  v=numeric(88*3)
  M=matrix(v,ncol=3)
  bx=c(lx,ux) 
  by=c(ly,uy)
  bxy=c(lxy,uxy)
  s1=function(bx,by,bxy) c(((-2*p+sqrt((2*p)^2-4*(-by*bxy)^2))/(-2*by*bxy))^2,by^2,bxy)
  s2=function(bx,by,bxy) c(((-2*p-sqrt((2*p)^2-4*(-by*bxy)^2))/(-2*by*bxy))^2,by^2,bxy)
  s3=function(bx,by,bxy) c((p+1)/(bxy+1),(p+1)/(bxy+1),bxy)
  s4=function(bx,by,bxy) c((p-1)/(bxy-1),(p-1)/(bxy-1),bxy)
  s5=function(bx,by,bxy) c(bx^2,by^2,bxy)
  s6=function(bx,by,bxy) c(bx,by,(p+sqrt(1-bx^2)*sqrt(1-by^2))/(bx*by))
  s7=function(bx,by,bxy) c(bx,by,(p-sqrt(1-bx^2)*sqrt(1-by^2))/(bx*by))
  s8=function(bx,by,bxy) c(bx^2,((-(-2*bx*bxy*p)+sqrt((-2*bx*bxy*p)^2-4*(bx^2*bxy^2+1-bx^2)*(bx^2-1+p^2)))/(2*(bx^2*bxy^2+1-bx^2)))^2,bxy)
  s9=function(bx,by,bxy) c(bx^2,((-(-2*bx*bxy*p)-sqrt((-2*bx*bxy*p)^2-4*(bx^2*bxy^2+1-bx^2)*(bx^2-1+p^2)))/(2*(bx^2*bxy^2+1-bx^2)))^2,bxy)
  s10=function(bx,by,bxy) c(((-(-2*bx*bxy*p)+sqrt((-2*bx*bxy*p)^2-4*(bx^2*bxy^2+1-bx^2)*(bx^2-1+p^2)))/(2*(bx^2*bxy^2+1-bx^2)))^2,by^2,bxy)
  s11=function(bx,by,bxy) c(((-(-2*bx*bxy*p)-sqrt((-2*bx*bxy*p)^2-4*(bx^2*bxy^2+1-bx^2)*(bx^2-1+p^2)))/(2*(bx^2*bxy^2+1-bx^2)))^2,by^2,bxy)
  w=c(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11)
  for (i in 1:2) {
    for (j in 1:2) {
      for (k in 1:2) {
        for (l in 1:11) {
          r=((l-1)*8)+(4*(i-1)+2*(j-1)+1*(k-1)+1)
          M[r,]=w[[l]](bx[i],by[j],bxy[k])
        }}}}
  M
  feas=function(arg) 
  {
    arg[1]>=(lx2-tol) & 
      arg[1]<=(ux2+tol) &
      arg[2]>=(ly2-tol) &
      arg[2]<=(uy2+tol) &
      arg[3]>=(lxy-tol) &
      arg[3]<=(uxy+tol) 
  }
  fs=apply(M,1,feas)
  N=M[which(fs==TRUE),]
  obj=function(arg) sr*(p-sqrt(arg[1])*sqrt(arg[2])*arg[3])/(1-arg[1])
  l=min(apply(N,1,obj), na.rm = TRUE)
  u=max(apply(N,1,obj), na.rm = TRUE)
  return(c(l,u))
}

BB.confound <- function(x,y,w) {
  start <- Sys.time()
  w <- data.frame(w)
  setClass("Node", slots = list(I = "vector", J = "ANY", beta = "ANY", depth = "numeric", x.res = "ANY", y.res = "ANY", J.res = "ANY"))
  
  max.depth <- ncol(w)
  nodes.visited <- rep(0,2)
  unique.nodes <- rep(0,2)
  confound.int <- NULL
  
  for (i in 1:2) {
    print(paste("Starting", c("min","max")[i]))
    node.stack <- list(new("Node", I = numeric(0), J = 1:max.depth, depth = 0, x.res = lm(x~1)$residuals, y.res = lm(y~1)$residuals, J.res = calc.J.res(w, numeric(0), 1:max.depth)))
    current.best <- node.stack[[1]]@beta <- lm(node.stack[[1]]@y.res~node.stack[[1]]@x.res)$coef[2] * (-1)^(i+1)
    while (length(node.stack) > 0) {
      nodes.visited[i] <- nodes.visited[i] + 1
      if (nodes.visited[i] %% 10000 == 0) {
        print(paste("Progress:", nodes.visited[i], "visited"))
      }
      # Pop first node off stack
      current <- node.stack[[1]]
      node.stack <- node.stack[-1]
      
      # Calculate beta if applicable
      if (is.null(current@beta)) {
        current@beta <- as.numeric(lm(current@y.res ~ current@x.res)$coef[2]) * (-1)^(i+1)
        current.best <- min(current.best, current@beta)
      }
      
      # Check to see if max.depth reached
      if (current@depth == max.depth) {
        next
      }
      
      # Calculate confounding interval and add new nodes if applicable
      # Save x.res, y.res, and J.res (minus one column) for right child since it is the same
      for (node in list(new("Node", I = c(current@I,current@J[1]), J = current@J[-1], depth = current@depth+1),
                        new("Node", I = current@I, J = current@J[-1], beta = current@beta, depth = current@depth+1, x.res = current@x.res, y.res = current@y.res, J.res = if(current@depth < max.depth - 1) current@J.res[,-1, drop = FALSE]))) {
        # Calculate x.res, y.res, and J.res for left child
        if (is.null(node@beta)) {
          unique.nodes[i] <- unique.nodes[i] + 1
          node@x.res <- lm(x~., data = w[node@I])$residuals # x ~ I
          node@y.res <- lm(y~., data = w[node@I])$residuals # y ~ I
          node@J.res <- calc.J.res(w, node@I, node@J) # J ~ I
        }
        
        # Calculate upper.r2.x, upper.r2.y
        if (!is.null(node@J.res)) {
          upper.r2.x <- summary(lm(node@x.res~., data = node@J.res))$r.squared
          upper.r2.y <- summary(lm(node@y.res~., data = node@J.res))$r.squared
        } else {
          upper.r2.x <- 0
          upper.r2.y <- 0
        }
        
        # Confounding Interval
        b <- suppressWarnings(f(cor(node@x.res,node@y.res), sd(node@y.res) / sd(node@x.res), 0, upper.r2.x, 0, upper.r2.y, -1, 1))[i] * (-1)^(i+1)
        if (b < current.best) {
          # Add node to stack if the lower bound is less than current.best
          node.stack <- c(node.stack, node)
        }
      }
    }
    confound.int <- c(confound.int, current.best * (-1)^(i+1))
  }
  end <- Sys.time()
  
  print(paste0("Find min: ", nodes.visited[1], " of ", (2^(max.depth+1)-1), " nodes visited (", round(nodes.visited[1] / (2^(max.depth+1)-1) * 100,2), "%)"))
  print(paste0("          ", unique.nodes[1], " of ", 2^max.depth, " unique nodes explored (",round(unique.nodes[1] / (2^max.depth) * 100,2), "%)"))
  print(paste0("Find max: ", nodes.visited[2], " of ", (2^(max.depth+1)-1), " nodes visited (", round(nodes.visited[2] / (2^(max.depth+1)-1) * 100,2), "%)"))
  print(paste0("          ", unique.nodes[2], " of ", 2^max.depth, " unique nodes explored (",round(unique.nodes[2] / (2^max.depth) * 100,2), "%)"))
  print(paste("Total Time elapsed:", as.numeric(difftime(end, start, units = "secs")), "seconds"))
  return(confound.int)
}


BB.confound.reorder <- function(x,y,w) {
  start <- Sys.time()
  w <- data.frame(w)
  setClass("Node", slots = list(I = "vector", J = "ANY", beta = "ANY", depth = "numeric", x.res = "ANY", y.res = "ANY", J.res = "ANY"))
  
  max.depth <- ncol(w)
  nodes.visited <- rep(0,2)
  unique.nodes <- rep(0,2)
  confound.int <- NULL
  
  for (i in 1:2) {
    print(paste("Starting", c("min","max")[i]))
    root <- new("Node", I = numeric(0), depth = 0, x.res = lm(x~1)$residuals, y.res = lm(y~1)$residuals, J.res = calc.J.res(w, numeric(0), 1:max.depth), J = 1:max.depth)
    root@J <- root@J[order(abs(cor(root@x.res, root@J.res) * cor(root@y.res, root@J.res)), decreasing = TRUE)]
    current.best <- root@beta <- lm(root@y.res~root@x.res)$coef[2] * (-1)^(i+1)
    node.stack <- list(root)
    while (length(node.stack) > 0) {
      nodes.visited[i] <- nodes.visited[i] + 1
      if (nodes.visited[i] %% 10000 == 0) {
        print(paste("Progress:", nodes.visited[i], "visited"))
      }
      # Pop first node off stack
      current <- node.stack[[1]]
      node.stack <- node.stack[-1]
      
      # Calculate beta.x and reorder variables if applicable
      if (is.null(current@beta)) {
        current@beta <- as.numeric(lm(current@y.res ~ current@x.res)$coef[2]) * (-1)^(i+1)
        current.best <- min(current.best, current@beta)
        if (length(current@J) > 1) {
          # Reorder variables for left child
          product <- abs(cor(current@x.res, current@J.res) * cor(current@y.res, current@J.res))
          next.J <- which(product == max(product))
          current@J <- c(current@J[next.J], current@J[-next.J])
          # current@J <- current@J[order(abs(cor(current@x.res, current@J.res) * cor(current@y.res, current@J.res)), decreasing = TRUE)]
        }
      }
      
      # Check to see if max.depth reached
      if (current@depth == max.depth) {
        next
      }
      
      # Calculate confounding interval and add new nodes if applicable
      # Save x.res, y.res, and J.res (minus one column) for right child since it is the same
      for (node in list(new("Node", I = c(current@I,current@J[1]), J = current@J[-1], depth = current@depth+1),
                        new("Node", I = current@I, J = current@J[-1], beta = current@beta, depth = current@depth+1, x.res = current@x.res, y.res = current@y.res, J.res = if(current@depth < max.depth - 1) current@J.res[,-1, drop = FALSE]))) {
        # Calculate x.res, y.res, and J.res for left child
        if (is.null(node@beta)) {
          unique.nodes[i] <- unique.nodes[i] + 1
          node@x.res <- lm(x~., data = w[node@I])$residuals # x ~ I
          node@y.res <- lm(y~., data = w[node@I])$residuals # y ~ I
          node@J.res <- calc.J.res(w, node@I, node@J) # J ~ I
        }
        
        # Calculate upper.r2.x, upper.r2.y
        if (!is.null(node@J.res)) {
          upper.r2.x <- summary(lm(node@x.res~., data = node@J.res))$r.squared
          upper.r2.y <- summary(lm(node@y.res~., data = node@J.res))$r.squared
        } else {
          upper.r2.x <- 0
          upper.r2.y <- 0
        }
        
        # Confounding Interval
        b <- suppressWarnings(f(cor(node@x.res,node@y.res), sd(node@y.res) / sd(node@x.res), 0, upper.r2.x, 0, upper.r2.y, -1, 1))[i] * (-1)^(i+1)
        if (b < current.best) {
          # Add node to stack if the lower bound is less than current.best
          node.stack <- c(node.stack, node)
        }
      }
    }
    confound.int <- c(confound.int, current.best * (-1)^(i+1))
  }
  end <- Sys.time()
  print(paste0("Find min: ", nodes.visited[1], " of ", (2^(max.depth+1)-1), " nodes visited (", round(nodes.visited[1] / (2^(max.depth+1)-1) * 100,2), "%)"))
  print(paste0("          ", unique.nodes[1], " of ", 2^max.depth, " unique nodes visited (",round(unique.nodes[1] / (2^max.depth) * 100,2), "%)"))
  print(paste0("Find max: ", nodes.visited[2], " of ", (2^(max.depth+1)-1), " nodes visited (", round(nodes.visited[2] / (2^(max.depth+1)-1) * 100,2), "%)"))
  print(paste0("          ", unique.nodes[2], " of ", 2^max.depth, " unique nodes visited (",round(unique.nodes[2] / (2^max.depth) * 100,2), "%)"))
  print(paste("Total Time elapsed:", as.numeric(difftime(end, start, units = "secs")), "seconds"))
  return(confound.int)
}