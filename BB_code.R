# This file provides open source software to support the methodologies of the paper
# Branch and Bound to Assess Stability of Regression Coefficients in Uncertain Models
# by Knaeble, Razo, Hughes, and Abramson

# See README for details

# Helper function to calculate residuals of I
calc.I.res <- function(w, J, I) {
  I.res <- matrix(nrow = nrow(w), ncol = length(I))
  if (length(J) == 0) { # Check for root node
    # Create matrix of residuals from centering vectors of I
    for (i in I) {
      I.res[,which(I == i)] <- lm(w[,i]~1)$residuals
    }
  } else if (length(I) == 0) { # Check for leaf node
    return(NULL)
  } else { # Intermediate node
    # Create matrix of residuals from regressing vectors of I onto span of J
    for (i in I) {
      I.res[,which(I == i)] <- lm(w[,i]~., data = w[J])$residuals
    }
  }
  return(data.frame(I.res))
}

# Confounding interval function supplied by Knaeble
# Found at https://github.com/bknaeble/ConfoundingIntervals/tree/master
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
  l=min(apply(N,1,obj), na.rm = TRUE) # I added "na.rm = TRUE"
  u=max(apply(N,1,obj), na.rm = TRUE) # I added "na.rm = TRUE"
  return(c(l,u))
}

# Branch and Bound (BB) Algorithm
BB.confound <- function(x,y,w) {
  start <- Sys.time()
  w <- data.frame(w)
  # Initialize "Node" class
  setClass("Node", slots = list(J = "ANY", I = "ANY", beta = "ANY", depth = "numeric",
                                x.res = "ANY", y.res = "ANY", I.res = "ANY"))
  # Initialize variables
  max.depth <- ncol(w)
  nodes.visited <- rep(0,2)
  confound.int <- NULL
  min.J <- NULL
  max.J <- NULL
  
  # iter_no = 1 for min, iter_no = 2 for max
  for (iter_no in 1:2) {
    print(paste0("Starting ", c("min","max")[iter_no], " (", Sys.time(), ")"))
    # Root node
    root <- new("Node", J = numeric(0), I = 1:max.depth, depth = 0, x.res = lm(x~1)$residuals,
                y.res = lm(y~1)$residuals, I.res = calc.I.res(w, numeric(0), 1:max.depth))
    current.best <- root@beta <- lm(root@y.res~root@x.res)$coef[2] * (-1)^(iter_no+1)
    queue <- list(root)
    while (length(queue) > 0) {
      nodes.visited[iter_no] <- nodes.visited[iter_no] + 1
      # Progress update
      if (nodes.visited[iter_no] %% 100 == 0) {
        print(paste0("Progress (", Sys.time(), "): ", nodes.visited[iter_no], " visited. ", length(queue), " on stack."))
      }
      
      # Take first node out of queue
      current <- queue[[1]]
      queue <- queue[-1]
      
      # Check to see if max.depth reached (leaf node)
      if (current@depth == max.depth) {
        next
      }
      
      # Check to see if we should add child nodes to the queue
      # Save x.res, y.res, and I.res (minus one column) for left child since it is the same
      for (node in list(new("Node", J = c(current@J,current@I[1]), I = current@I[-1], depth = current@depth+1),
                        new("Node", J = current@J, I = current@I[-1], beta = current@beta, depth = current@depth+1,
                            x.res = current@x.res, y.res = current@y.res,
                            I.res = if(current@depth < max.depth - 1) current@I.res[,-1, drop = FALSE]))) {
        # Calculate x.res, y.res, and I.res for right child
        if (is.null(node@beta)) {
          node@x.res <- lm(x~., data = w[node@J])$residuals
          node@y.res <- lm(y~., data = w[node@J])$residuals
          node@I.res <- calc.I.res(w, node@J, node@I)
        }
        
        # Calculate r2.wI.x, r2.wI.y
        if (!is.null(node@I.res)) {
          r2.wI.x <- summary(lm(node@x.res~., data = node@I.res))$r.squared
          r2.wI.y <- summary(lm(node@y.res~., data = node@I.res))$r.squared
        } else {
          r2.wI.x <- 0
          r2.wI.y <- 0
        }
        
        # Confounding Interval
        b <- suppressWarnings(f(cor(node@x.res,node@y.res), sd(node@y.res) / sd(node@x.res),
                                0, r2.wI.x, 0, r2.wI.y, -1, 1))[iter_no] * (-1)^(iter_no+1)
        if (b < current.best) {
          # Calculate beta.x and update current.best if applicable
          if (is.null(node@beta)) {
            node@beta <- as.numeric(lm(node@y.res ~ node@x.res)$coef[2]) * (-1)^(iter_no+1)
            if (current.best > node@beta) {
              if (iter_no == 1) {
                min.J <- node@J
              } else {
                max.J <- node@J
              }
            }
            current.best <- min(current.best, node@beta)
          }
          # Put node in queue if the lower bound is less than current.best
          queue <- c(queue, node)
        }
      }
    }
    confound.int <- c(confound.int, current.best * (-1)^(iter_no+1))
  }
  end <- Sys.time()
  print("------------------------")
  print(paste0("Find min: ", nodes.visited[1], " of ", (2^(max.depth+1)-1), " nodes visited (", round(nodes.visited[1] / (2^(max.depth+1)-1) * 100,2), "%)"))
  print(paste("Min vars:", paste(colnames(w)[sort(min.J)], collapse = ",")))
  print("----------")
  print(paste0("Find max: ", nodes.visited[2], " of ", (2^(max.depth+1)-1), " nodes visited (", round(nodes.visited[2] / (2^(max.depth+1)-1) * 100,2), "%)"))
  print(paste("Max vars:", paste(colnames(w)[sort(max.J)], collapse = ",")))
  print("----------")
  print(paste(sum(nodes.visited), "total nodes visited."))
  print(paste("Total Time elapsed:", as.numeric(difftime(end, start, units = "secs")), "seconds"))
  print("------------------------")
  return(confound.int)
}

# Branch and Bound with Reordering (BBR) Algorithm
BB.confound.reorder <- function(x,y,w) {
  start <- Sys.time()
  w <- data.frame(w)
  # Initialize "Node" class
  setClass("Node", slots = list(J = "ANY", I = "ANY", beta = "ANY", depth = "numeric",
                                x.res = "ANY", y.res = "ANY", I.res = "ANY"))  
  # Initialize variables
  max.depth <- ncol(w)
  nodes.visited <- rep(0,2)
  confound.int <- NULL
  min.J <- NULL
  max.J <- NULL
  
  # iter_no = 1 for min, iter_no = 2 for max
  for (iter_no in 1:2) {
    print(paste0("Starting ", c("min","max")[iter_no], " (", Sys.time(), ")"))
    # Root node
    root <- new("Node", J = numeric(0), I = 1:max.depth, depth = 0, x.res = lm(x~1)$residuals,
                y.res = lm(y~1)$residuals, I.res = calc.I.res(w, numeric(0), 1:max.depth))
    root@I <- root@I[order(abs(cor(root@x.res, root@I.res) * cor(root@y.res, root@I.res)), decreasing = TRUE)]
    current.best <- root@beta <- lm(root@y.res~root@x.res)$coef[2] * (-1)^(iter_no+1)
    queue <- list(root)
    while (length(queue) > 0) {
      nodes.visited[iter_no] <- nodes.visited[iter_no] + 1
      # Progress update
      if (nodes.visited[iter_no] %% 100 == 0) {
        print(paste0("Progress (", Sys.time(), "): ", nodes.visited[iter_no], " visited. ", length(queue), " on stack."))
      }
      
      # Take first node out of queue
      current <- queue[[1]]
      queue <- queue[-1]
      
      # Check to see if max.depth reached (leaf node)
      if (current@depth == max.depth) {
        next
      }
      
      # Check to see if we should add child nodes to the queue
      # Save x.res, y.res, and I.res (minus one column) for left child since it is the same
      for (node in list(new("Node", J = c(current@J,current@I[1]), I = current@I[-1], depth = current@depth+1),
                        new("Node", J = current@J, I = current@I[-1], beta = current@beta, depth = current@depth+1,
                            x.res = current@x.res, y.res = current@y.res,
                            I.res = if(current@depth < max.depth - 1) current@I.res[,-1, drop = FALSE]))) {
        # Calculate x.res, y.res, and I.res for right child
        if (is.null(node@beta)) {
          node@x.res <- lm(x~., data = w[node@J])$residuals
          node@y.res <- lm(y~., data = w[node@J])$residuals
          node@I.res <- calc.I.res(w, node@J, node@I)
        }
        
        # Calculate r2.wI.x, r2.wI.y
        if (!is.null(node@I.res)) {
          r2.wI.x <- summary(lm(node@x.res~., data = node@I.res))$r.squared
          r2.wI.y <- summary(lm(node@y.res~., data = node@I.res))$r.squared
        } else {
          r2.wI.x <- 0
          r2.wI.y <- 0
        }
        
        # Confounding Interval
        b <- suppressWarnings(f(cor(node@x.res,node@y.res), sd(node@y.res) / sd(node@x.res),
                                0, r2.wI.x, 0, r2.wI.y, -1, 1))[iter_no] * (-1)^(iter_no+1)
        if (b < current.best) {
          # Calculate beta.x and update current.best if applicable
          if (is.null(node@beta)) {
            node@beta <- as.numeric(lm(node@y.res ~ node@x.res)$coef[2]) * (-1)^(iter_no+1)
            if (current.best > node@beta) {
              if (iter_no == 1) {
                min.J <- node@J
              } else {
                max.J <- node@J
              }
            }
            current.best <- min(current.best, node@beta)
          }
          
          # Reorder I covariates
          if (length(node@I) > 1) {
            product <- abs(cor(node@x.res, node@I.res) * cor(node@y.res, node@I.res))
            next.I <- which(product == max(product))
            node@I <- c(node@I[next.I], node@I[-next.I])
          }
          
          # Put node in queue if the lower bound is less than current.best
          queue <- c(queue, node)
        }
      }
    }
    confound.int <- c(confound.int, current.best * (-1)^(iter_no+1))
  }
  end <- Sys.time()
  print("------------------------")
  print(paste0("Find min: ", nodes.visited[1], " of ", (2^(max.depth+1)-1), " nodes visited (", round(nodes.visited[1] / (2^(max.depth+1)-1) * 100,2), "%)"))
  print(paste("Min vars:", paste(colnames(w)[sort(min.J)], collapse = ",")))
  print("----------")
  print(paste0("Find max: ", nodes.visited[2], " of ", (2^(max.depth+1)-1), " nodes visited (", round(nodes.visited[2] / (2^(max.depth+1)-1) * 100,2), "%)"))
  print(paste("Max vars:", paste(colnames(w)[sort(max.J)], collapse = ",")))
  print("----------")
  print(paste(sum(nodes.visited), "total nodes visited."))
  print(paste("Total Time elapsed:", as.numeric(difftime(end, start, units = "secs")), "seconds"))
  print("------------------------")
  return(confound.int)
}

